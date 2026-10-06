/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file Voronoi.h
 * @brief Spherical Voronoi: animated sites shaded by nearest-site distance,
 *        with edge sharpening and optional cell borders.
 */

#include "core/engine/engine.h"
#include "core/spatial/kd_tree.h"

#include <cmath>
#include <span>

// Unit-test accessor for the seeded sites and the coarse-grid tuning constants.
namespace hs_test {
namespace effects_tests {
struct VoronoiWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Spherical Voronoi effect: scatters sites on the unit sphere, animates
 *        them, and shades each pixel by its nearest site with edge-sharpening
 *        and optional cell borders.
 * @tparam W Framebuffer width in pixels.
 * @tparam H Framebuffer height in pixels.
 */
template <int W, int H> class Voronoi : public Effect {
public:
  static constexpr const char *EFFECT_ID = "Voronoi";

  /** @brief Construction config. */
  static constexpr EffectConfig CONFIG{.strobe = true};

  /**
   * @brief Constructs the effect with the templated framebuffer dimensions.
   */
  HS_COLD_MEMBER Voronoi() : Effect(W, H, CONFIG) {}

  /**
   * @brief Configures arenas, registers GUI params, allocates the sites buffer,
   *        and seeds the initial sites.
   */
  HS_COLD_MEMBER void init() override {
    // Persistent holds the sites buffer; scratch_arena_a holds the per-frame
    // KD-tree (positions + nodes + build indices).
    ArenaSplit{SCRATCH_A_BYTES, 0}.configure();

    register_int_param("Num Sites", &params.num_sites, 1, MAX_SITES);
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
    Canvas canvas(*this);

    {
      HS_PROFILE(vo_animate);
      if (active_site_count() != current_num_sites)
        seed_sites();

      float s = logf(params.speed + 1.0f) * SITE_SPIN_RADIANS;
      // Every site turns by the same angle, so the half-angle trig is shared;
      // only the axis differs. Expands make_rotation(site.axis, s) exactly.
      const float half_cos = cosf(s / 2);
      const float half_sin = sinf(s / 2);

      for (size_t i = 0; i < sites_buffer.size(); ++i) {
        auto &site = sites_buffer[i];
        // Renormalize: rotate() drifts |pos| off the unit sphere over a long run,
        // and the unit-site invariants (nearest-by-Euclidean == nearest-by-max-dot,
        // and the border acosf(dot) staying in range) require unit vectors.
        const math::Quaternion q(half_cos, half_sin * site.axis);
        site.pos = math::rotate(site.pos, q).normalized();
      }
    }

    // Build a KD-tree over the moving site positions once per frame. On the unit
    // sphere nearest-by-Euclidean == nearest-by-max-dot (|p-s|^2 = 2 - 2*p*s),
    // so the k=2 query is exact.
    ScratchScope scope_guard(scratch_arena_a);
    math::Vector *positions =
        scratch_arena_a.allocate_n<math::Vector>(sites_buffer.size());
    // IIFE times the per-frame KD build while keeping `tree` at frame scope
    // (guaranteed copy elision on the prvalue return).
    KDTree tree = [&]() -> KDTree {
      HS_PROFILE(vo_kdtree);
      for (size_t i = 0; i < sites_buffer.size(); ++i)
        positions[i] = sites_buffer[i].pos;
      return KDTree(scratch_arena_a, std::span<const math::Vector>(
                                         positions, sites_buffer.size()));
    }();

    // One node for all per-pixel work (corner pre-pass + shading loop); never
    // scope inside the per-pixel loop.
    HS_PROFILE(vo_shade);

    // Coarse-grid coherence: classify the nearest pair once per coarse-grid
    // corner, then shade every pixel of a block from the deduped union of its
    // four corners' pairs (<= 8 candidate sites) by an exact top-2 dot scan.
    // A cell missed by all four corners is dropped.
    // Voronoi cell pixel size falls as ~1/sqrt(num_sites), so shrink the block
    // with the site count, floored at the edge MAX_SITES would give at this H.
    // Full-sphere row-pitch estimate: an equal-area cell spans sqrt(4π/n)·H/π
    // rows, and this edge is 1/sqrt(π) ≈ 0.56 of that estimate.
    const float cell_px = (2.0f * H / math::PI_F) /
                          sqrtf(static_cast<float>(sites_buffer.size()));
    const int B = hs::clamp(static_cast<int>(cell_px), COHERENCE_BLOCK_MIN,
                            COHERENCE_BLOCK);
    if (B == 1) {
      Scan::Shader::draw<W, H>(canvas, [&](const math::Vector &p) {
        const auto nearest = tree.nearest(p, 2);
        const uint16_t i0 = nearest[0].original_index;
        const uint16_t i1 = nearest.size() > 1 ? nearest[1].original_index : i0;
        return shade(i0, math::dot(p, sites_buffer[i0].pos), i1,
                     nearest.size() > 1 ? math::dot(p, sites_buffer[i1].pos)
                                        : NO_DOT);
      });
      return;
    }
    Scan::Shader::draw_block_coherent<W, H, SITES_PER_CORNER>(
        canvas, B, positions, scratch_arena_a,
        [&](const math::Vector &p) {
          const CellId pair = classify(tree, p);
          return std::array<uint16_t, SITES_PER_CORNER>{pair.lo, pair.hi};
        },
        [&](const math::Vector &p, const CandSet &cs)
            __attribute__((always_inline)) {
              float d0 = NO_DOT, d1 = NO_DOT;
              uint8_t b0 = 0, b1 = 0;
              for (uint8_t i = 0; i < cs.n; ++i) {
                float d = math::dot(p, cs.pos[i]);
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
              return shade(cs.idx[b0], d0, cs.idx[b1], d1);
            });
  }

private:
  // Test seam: reaches the seeded sites and the coarse-grid tuning constants.
  friend struct ::hs_test::effects_tests::VoronoiWhiteBox;

  /**
   * @brief One Voronoi seed: position on the unit sphere, the rotation axis it
   *        spins about, and the cell fill color.
   */
  struct Site {
    math::Vector pos;  /**< Position on the unit sphere. */
    math::Vector axis; /**< Unit rotation axis the site spins about. */
    Color4 color;      /**< Cell fill color (16-bit linear channels). */
  };

  static constexpr int MAX_SITES = 400; /**< Buffer capacity; the sites buffer
                                             is allocated once at this size. */
  static constexpr float SITE_SPIN_RADIANS =
      0.005f; /**< Radians per frame after the logarithmic speed mapping. */
  static constexpr int COHERENCE_BLOCK = 8; /**< Coarse-coherence block edge in
      pixels at low site counts: each pixel shades from the union of its block
      corners' nearest pairs. Smaller is safer (fewer missed sub-block cells)
      but classifies more corners; the render path shrinks the block toward
      COHERENCE_BLOCK_MIN as the site count rises. */
  static constexpr int COHERENCE_BLOCK_MIN = std::max(
      1, static_cast<int>(
             (2.0 * H / static_cast<double>(math::PI_F)) /
             math::constexpr_sqrt(MAX_SITES))); /**< Smallest adaptive block
      edge: the render path's edge formula evaluated at MAX_SITES, floored at
      one pixel. */
  static_assert(COHERENCE_BLOCK_MIN <= COHERENCE_BLOCK,
                "Voronoi coherence block bounds are inverted");

  /** @brief Canonical nearest-pair identity at a sample point. */
  struct CellId {
    uint16_t lo; /**< min(nearest, second) site index. */
    uint16_t hi; /**< max(nearest, second) site index. */
  };

  /** @brief Below any unit-vector dot product, so it doubles as the "no
      candidate" marker the top-2 scan hands to shade(). */
  static constexpr float NO_DOT = -2.0f;

  /**
   * @brief Canonical (order-independent) nearest-pair identity at a sample
   *        point.
   * @param tree KD-tree over the current frame's site positions.
   * @param p Sample point on the unit sphere.
   * @return The nearest pair as an ordered {lo, hi} index set, so the two
   *         query orders along a cell seam map to the same identity.
   */
  static CellId classify(const KDTree &tree, const math::Vector &p) {
    auto knn = tree.nearest(p, 2);
    uint16_t a = knn[0].original_index;
    uint16_t b = knn.size() > 1 ? knn[1].original_index : a;
    return {std::min(a, b), std::max(a, b)};
  }

  /**
   * @brief Resolves one pixel's color from its already-identified nearest
   *        pair.
   * @param i0 Nearest site index.
   * @param d0 Dot with the nearest site, so the larger of the two.
   * @param i1 Second-nearest site index.
   * @param d1 Dot with the second site; NO_DOT when the block held a single
   *        candidate.
   * @return The cell color, alpha 0 on a border seam.
   */
  __attribute__((always_inline)) Color4 shade(uint16_t i0, float d0,
                                              uint16_t i1, float d1) const {
    const Site &best_site = sites_buffer[i0];
    const bool has_second = d1 > NO_DOT;

    Color4 c = best_site.color;

    // Border sharpening: a larger sharpness saturates `factor` for smaller
    // nearest/second-nearest gaps, shrinking the cross-cell blend band.
    if (has_second && params.sharpness > 0.0f) {
      const Site &sec_site = sites_buffer[i1];
      float diff = d0 - d1;
      float factor = fminf(1.0f, diff * params.sharpness);
      factor = math::quintic_kernel(factor);
      float t = 0.5f + 0.5f * factor;

      uint16_t frac = static_cast<uint16_t>(t * 65535.0f + 0.5f);
      c.color = sec_site.color.color.lerp16(best_site.color.color, frac);
    }

    // A "Border Thick" of 0 skips the two acosf calls. d0 is the nearest, so
    // d0 >= d1 and the cell gap is non-negative.
    if (params.border_thickness > 0.0f && has_second) {
      float dist1 = acosf(hs::clamp(d0, -1.0f, 1.0f));
      float dist2 = acosf(hs::clamp(d1, -1.0f, 1.0f));
      if (dist2 - dist1 < params.border_thickness) {
        // Paint the seam black. The per-pixel store writes color*alpha, so an
        // alpha-0 sample collapses to (0,0,0).
        c = Color4(0, 0, 0, 0);
      }
    }

    return c;
  }

  static constexpr int SITES_PER_CORNER = 2;
  using CandSet = Scan::Shader::BlockCandidates<SITES_PER_CORNER>;

  // Compile-time high-water check for the 64 KB scratch_arena_a reserve. Two
  // transient peaks share that arena, with positions + KD nodes live across both:
  //   build:   positions + KD nodes + KD build-index scratch
  //   shading: positions + KD nodes + corner cells + candidate-row sets
  static constexpr size_t SCRATCH_A_BYTES = 64 * 1024;
  static constexpr size_t POSITIONS_BYTES =
      size_t(MAX_SITES) * sizeof(math::Vector);
  static constexpr size_t KD_NODES_BYTES = size_t(MAX_SITES) * sizeof(KDNode);
  static constexpr size_t KD_BUILD_SCRATCH_BYTES =
      size_t(MAX_SITES) * sizeof(int);
  static constexpr size_t WALKER_BYTES =
      Scan::Shader::block_coherent_scratch_bytes<W, H, SITES_PER_CORNER,
                                                 COHERENCE_BLOCK_MIN>();
  static constexpr size_t SCRATCH_HIGH_WATER =
      POSITIONS_BYTES + KD_NODES_BYTES +
      (KD_BUILD_SCRATCH_BYTES > WALKER_BYTES ? KD_BUILD_SCRATCH_BYTES
                                             : WALKER_BYTES);
  static_assert(
      SCRATCH_HIGH_WATER <= SCRATCH_A_BYTES,
      "Voronoi scratch_arena_a budget too small for MAX_SITES positions + "
      "KD-tree + coarse-grid cells; raise SCRATCH_A_BYTES or lower MAX_SITES / "
      "coarsen COHERENCE_BLOCK");

  // init() binds the MAX_SITES sites buffer from the persistent arena, which
  // configure_arenas() sizes as the global arena less SCRATCH_A_BYTES.
  static constexpr size_t FOOTPRINT_BYTES =
      size_t(MAX_SITES) * sizeof(Site) + alignof(Site);
  static_assert(FOOTPRINT_BYTES <=
                    ArenaSplit{SCRATCH_A_BYTES, 0}.device_persistent(),
                "Voronoi persistent footprint exceeds its device partition; "
                "lower MAX_SITES or shrink SCRATCH_A_BYTES");

  int current_num_sites = 0;      /**< Count currently seeded; re-seeds (clear +
                                   refill, no realloc) when the slider changes. */
  ArenaVector<Site> sites_buffer; /**< Active Voronoi sites for the frame. */

  /**
   * @brief Returns the active site count from the "Num Sites" slider.
   * @return Integer slider value clamped to [1, MAX_SITES].
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

    for (int i = 0; i < n; i++) {
      math::Vector pos = math::fib_spiral(n, /*eps=*/0.5f, i);

      math::Vector axis = math::random_vector();

      float t = i / (float)(n > 1 ? n - 1 : 1);
      Color4 color = Palettes::RICH_SUNSET.get(t);

      sites_buffer.push_back({pos, axis, color});
    }
    current_num_sites = n;
  }

  /**
   * @brief Live-tunable GUI parameters for the Voronoi effect.
   */
  struct Params {
    int num_sites = 200;      /**< Live-tunable site count (GUI slider). */
    float speed = 20.0f;      /**< Site spin rate (GUI slider). */
    float sharpness = 100.0f; /**< Edge sharpening; larger narrows the border
                                  blend band; 0 gives hard edges, like infinite sharpness. */
    float border_thickness = 0.0f; /**< Cell-seam border width; 0 disables. */
  } params;
};
