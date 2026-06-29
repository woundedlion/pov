/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"
#include "core/spatial.h"

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
  FLASHMEM Voronoi() : Effect(W, H) {}

  /**
   * @brief Configures arenas, registers GUI params, allocates the sites buffer,
   *        and seeds the initial sites.
   */
  void init() override {
    // Persistent holds the sites buffer; scratch_arena_a holds the per-frame
    // KD-tree (positions + nodes + build indices).
    configure_arenas(GLOBAL_ARENA_SIZE - 64 * 1024, 64 * 1024, 0);

    registerParam("Num Sites", &params.num_sites, 1.0f,
                  static_cast<float>(MAX_SITES));
    registerParam("Speed", &params.speed, 0.0f, 100.0f);
    registerParam("Sharpness", &params.sharpness, 1.0f, 500.0f);
    registerParam("Border Thick", &params.borderThickness, 0.0f, 0.1f);

    sites_buffer.bind(persistent_arena, MAX_SITES);
    seed_sites();
  }

  /// POV column-strobe flag; strobes (see Effect::strobe_columns).
  bool strobe_columns() const override { return true; }

  /**
   * @brief Animates the sites, builds a per-frame KD-tree, and shades each
   *        pixel by its nearest site (with edge sharpening and optional
   *        borders).
   */
  void draw_frame() override {
    Canvas canvas(*this);

    // Re-seed when the GUI changes the site count (integer change only, so
    // dragging within a bucket doesn't thrash).
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

    // Build a KD-tree over the moving site positions once per frame. On the unit
    // sphere nearest-by-Euclidean == nearest-by-max-dot (|p-s|^2 = 2 - 2*p*s),
    // so the k=2 query is exact.
    ScratchScope scope_guard(scratch_arena_a);
    Vector *positions = static_cast<Vector *>(scratch_arena_a.allocate(
        sites_buffer.size() * sizeof(Vector), alignof(Vector)));
    for (size_t i = 0; i < sites_buffer.size(); ++i)
      positions[i] = sites_buffer[i].pos;
    KDTree tree(scratch_arena_a,
                std::span<const Vector>(positions, sites_buffer.size()));

    // Resolves the final color from the already-identified nearest pair: the
    // nearest site index i0 and its dot d0 (the larger), and — when a second
    // neighbor exists — the second site index i1 and its dot d1.
    auto shade = [&](int i0, float d0, bool hasSecond, int i1,
                     float d1) -> Color4 {
      const Site &bestSite = sites_buffer[i0];
      float maxDot1 = d0;
      float maxDot2 = hasSecond ? d1 : -2.0f;

      Color4 c = bestSite.color;

      // Border sharpening: a larger sharpness saturates `factor` for smaller
      // nearest/second-nearest gaps, shrinking the cross-cell blend band.
      if (hasSecond && params.sharpness > 0.0f) {
        const Site &secSite = sites_buffer[i1];
        float diff = maxDot1 - maxDot2;
        float factor = std::min(1.0f, diff * params.sharpness);
        factor = quintic_kernel(factor);
        float t = 0.5f + 0.5f * factor;

        uint16_t frac = static_cast<uint16_t>(t * 65535.0f + 0.5f);
        c.color = secSite.color.color.lerp16(bestSite.color.color, frac);
      }

      // Borders — driven entirely by the "Border Thick" slider: a thickness of
      // 0 disables them (and skips the two acosf calls below), any positive
      // value paints the seam between the nearest two sites. d0 is the nearest,
      // so maxDot1 >= maxDot2 → dist1 <= dist2 and the cell gap is non-negative.
      if (params.borderThickness > 0.0f && hasSecond) {
        float dist1 = acosf(hs::clamp(maxDot1, -1.0f, 1.0f));
        float dist2 = acosf(hs::clamp(maxDot2, -1.0f, 1.0f));
        if (dist2 - dist1 < params.borderThickness) {
          // Paint the seam black. The Scan sink writes color*alpha, so an
          // alpha-0 fragment collapses to (0,0,0) regardless of its RGB.
          c = Color4(0, 0, 0, 0);
        }
      }

      return c;
    };

    // Canonical (order-independent) nearest-pair identity at a sample point: the
    // two query orders along a cell seam (best/second swap) map to the same
    // {lo, hi} set, so corner comparisons below are seam-stable.
    auto classify = [&](const Vector &p) -> CellId {
      auto knn = tree.nearest(p, 2);
      uint16_t a = knn[0].original_index;
      bool hasSecond = knn.size() > 1;
      uint16_t b = hasSecond ? knn[1].original_index : a;
      return {std::min(a, b), std::max(a, b), hasSecond};
    };
    auto same_cell = [](const CellId &x, const CellId &y) {
      return x.lo == y.lo && x.hi == y.hi && x.hasSecond == y.hasSecond;
    };

    // Coarse-grid coherence: the nearest-pair identity is piecewise-constant
    // over each Voronoi cell. Classify the pair once per coarse-grid corner; an
    // interior pixel whose four corners agree skips the k=2 query and shades from
    // two exact dot products, while a pixel whose corners disagree (a seam, or
    // the dense region near a pole) falls back to the full query. A cell smaller
    // than the block missed by all four corners is dropped.
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
    const int B = hs::clamp(static_cast<int>(cell_px), kCoherenceBlockMin,
                            kCoherenceBlock);
    const int nbx = (x1 - x0 - 1) / B + 2; // corner columns spanning [x0, x1)
    const int nby = (y1 - y0 - 1) / B + 2; // corner rows spanning    [y0, y1)
    CellId *cells = static_cast<CellId *>(
        scratch_arena_a.allocate(nbx * nby * sizeof(CellId), alignof(CellId)));
    // Clamp to the last valid pixel (x1/y1 are exclusive) so every classified
    // point indexes the trig LUT.
    auto corner_x = [&](int j) { return std::min(x0 + j * B, x1 - 1); };
    auto corner_y = [&](int k) { return std::min(y0 + k * B, y1 - 1); };
    for (int k = 0; k < nby; ++k)
      for (int j = 0; j < nbx; ++j)
        cells[k * nbx + j] =
            classify(pixel_to_vector<W, H>(corner_x(j), corner_y(k)));

    // SAMPLES=1: one sample at pixel center — the per-pixel KD query is too
    // heavy to supersample. Mirrors Scan::Shader::draw<W, H, 1>'s clip iteration
    // and telemetry (open-coded only because the coarse grid needs the integer
    // pixel coordinates the generic shader callback does not expose).
    Scan::ScopedRenderTimer timer_guard(canvas);
    for (int y = y0; y < y1; ++y) {
      const int ky = (y - y0) / B;
      for (int x = x0; x < x1; ++x) {
        const int jx = (x - x0) / B;
        const CellId &c00 = cells[ky * nbx + jx];
        const CellId &c10 = cells[ky * nbx + (jx + 1)];
        const CellId &c01 = cells[(ky + 1) * nbx + jx];
        const CellId &c11 = cells[(ky + 1) * nbx + (jx + 1)];

        Vector p = pixel_to_vector<W, H>(x, y);
        Color4 sample;
        if (same_cell(c00, c10) && same_cell(c00, c01) &&
            same_cell(c00, c11)) {
          // Corners agree: reuse the pair identity, recompute the two exact
          // dots and order them (nearest = larger dot) — bit-identical to the
          // full query whenever no third site is actually nearest here.
          float da = dot(p, sites_buffer[c00.lo].pos);
          if (c00.hasSecond) {
            float db = dot(p, sites_buffer[c00.hi].pos);
            sample = (da >= db) ? shade(c00.lo, da, true, c00.hi, db)
                                : shade(c00.hi, db, true, c00.lo, da);
          } else {
            sample = shade(c00.lo, da, false, c00.lo, da);
          }
        } else {
          // Cell seam / dense block: full k=2 query.
          auto knn = tree.nearest(p, 2);
          int i0 = knn[0].original_index;
          float d0 = dot(p, knn[0].point);
          if (knn.size() > 1)
            sample = shade(i0, d0, true, knn[1].original_index,
                           dot(p, knn[1].point));
          else
            sample = shade(i0, d0, false, i0, d0);
        }
        canvas(x, y) = sample.color * sample.alpha;
      }
    }
  }

  /**
   * @brief Live-tunable GUI parameters for the Voronoi effect.
   */
  struct Params {
    float num_sites = 200.0f;     /**< Live-tunable site count (GUI slider). */
    float speed = 20.0f;          /**< Site spin rate (GUI slider). */
    float sharpness = 100.0f;     /**< Edge sharpening; larger narrows the
                                       border blend band. */
    float borderThickness = 0.0f; /**< Cell-seam border width; 0 disables. */
  } params;

  static constexpr int MAX_SITES = 400; /**< Buffer capacity; the sites buffer
                                             is allocated once at this size. */
  static constexpr int kCoherenceBlock = 8; /**< Coarse-coherence block edge in
      pixels at low site counts: a pixel whose surrounding block corners agree on
      the nearest pair skips the per-pixel KD query. Smaller is safer (fewer
      missed sub-block cells) but skips fewer queries; the render path shrinks the
      block toward kCoherenceBlockMin as the site count rises. */
  static constexpr int kCoherenceBlockMin = 4; /**< Smallest adaptive block edge.
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
    bool hasSecond; /**< Whether a second neighbor exists (>= 2 sites). */
  };

  // Compile-time high-water check for the 64 KB scratch_arena_a reserve. Two
  // transient peaks share that arena, with positions + KD nodes live across both:
  //   build:   positions + KD nodes + KD build-index scratch
  //   shading: positions + KD nodes + coarse-grid corner cells
  static constexpr size_t kScratchABytes = 64 * 1024;
  static constexpr size_t kPositionsBytes = size_t(MAX_SITES) * sizeof(Vector);
  static constexpr size_t kKdNodesBytes = size_t(MAX_SITES) * sizeof(KDNode);
  static constexpr size_t kKdBuildScratchBytes = size_t(MAX_SITES) * sizeof(int);
  // Corner grid spans the full canvas in block-px steps, +2 for the inclusive
  // [0,W)/[0,H) end corners (mirrors render()'s nbx/nby at full clip). Sized at
  // kCoherenceBlockMin — the worst case (most corners) the adaptive block hits.
  static constexpr size_t kCornerCols = size_t((W - 1) / kCoherenceBlockMin + 2);
  static constexpr size_t kCornerRows = size_t((H - 1) / kCoherenceBlockMin + 2);
  static constexpr size_t kCellsBytes = kCornerCols * kCornerRows * sizeof(CellId);
  static constexpr size_t kScratchHighWater =
      kPositionsBytes + kKdNodesBytes +
      (kKdBuildScratchBytes > kCellsBytes ? kKdBuildScratchBytes : kCellsBytes);
  static_assert(kScratchHighWater <= kScratchABytes,
                "Voronoi scratch_arena_a budget (64 KB) too small for MAX_SITES "
                "positions + KD-tree + coarse-grid cells; raise the reserve in "
                "init() or lower MAX_SITES / coarsen kCoherenceBlock");

  int current_num_sites = 0; /**< Count currently seeded; re-seeds (clear +
                                  refill, no realloc) when the slider changes. */
  ArenaVector<Site> sites_buffer; /**< Active Voronoi sites for the frame. */

  /**
   * @brief Returns the active site count from the "Num Sites" slider.
   * @return Slider value rounded to an integer and clamped to [1, MAX_SITES].
   */
  int active_site_count() const {
    return std::clamp(static_cast<int>(params.num_sites), 1, MAX_SITES);
  }

  /**
   * @brief (Re)seeds the active sites for the current "Num Sites" slider value.
   * @details Clears and refills up to the active count (no re-allocation). Sites
   *          are placed via the Fibonacci-sphere distribution for an even
   *          spread, each with a random spin axis and a palette color.
   */
  void seed_sites() {
    const int n = active_site_count();
    sites_buffer.clear();

    const float goldenAngle = PI_F * (3.0f - sqrtf(5.0f));

    for (int i = 0; i < n; i++) {
      // Guard n == 1: a 0 denominator would give NaN y. A single site sits at
      // the pole (y = 1).
      int span = n > 1 ? n - 1 : 1;
      float y = 1.0f - (i / (float)span) * 2.0f;
      float radius = sqrtf(std::max(0.0f, 1.0f - y * y));
      float theta = goldenAngle * i;

      float x = cosf(theta) * radius;
      float z = sinf(theta) * radius;

      Vector pos = Vector(x, y, z);

      float rx = hs::rand_f() * 2.0f - 1.0f;
      float ry = hs::rand_f() * 2.0f - 1.0f;
      float rz = hs::rand_f() * 2.0f - 1.0f;
      Vector axis = Vector(rx, ry, rz).normalized();

      float t = i / (float)(n > 1 ? n - 1 : 1);
      Color4 color = Palettes::richSunset.get(t);

      sites_buffer.push_back({pos, axis, color});
    }
    current_num_sites = n;
  }
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Voronoi)
