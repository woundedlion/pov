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

#include "core/platform/effect_identity.h"
#include "core/engine/engine.h"
#include "core/render/temporal.h"
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
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY =
      hs::effect_identities::Voronoi::VALUE;
  /// Canonical machine ID.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;

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
    previous_positions.bind(persistent_arena, MAX_SITES);
    corner_pairs = persistent_arena.allocate_n<CornerPair>(CORNER_CAPACITY);
    seed_sites();
  }

  /**
   * @brief Animates the sites, builds a per-frame KD-tree, and shades each
   *        pixel by its nearest site (with edge sharpening and optional
   *        borders).
   * @details With column-time sampling, each column shades the sites at its
   * sweep fraction of the frame's spin from their previous positions.
   */
  void draw_frame() override {
    Canvas canvas(*this);

    advance_sites();

    // Build a KD-tree over the moving site positions once per frame. On the unit
    // sphere nearest-by-Euclidean == nearest-by-max-dot (|p-s|^2 = 2 - 2*p*s),
    // so the k=2 query is exact.
    ScratchScope scope_guard(scratch_arena_a);
    math::Vector *positions =
        scratch_arena_a.allocate_n<math::Vector>(sites_buffer.size());
    // IIFE scopes the profile node to the KD build; `tree` lives at frame scope.
    KDTree tree = [&]() -> KDTree {
      HS_PROFILE(vo_kdtree);
      for (size_t i = 0; i < sites_buffer.size(); ++i)
        positions[i] = sites_buffer[i].pos;
      return KDTree(scratch_arena_a, std::span<const math::Vector>(
                                         positions, sites_buffer.size()));
    }();

    if (Temporal::temporal_sampling_enabled() && frame_spin != 0.0f) {
      draw_temporal(canvas, tree);
      return;
    }

    corners_valid = false;

    // One node for all per-pixel work (corner pre-pass + shading loop); never
    // scope inside the per-pixel loop.
    HS_PROFILE(vo_shade);

    // Coarse-grid coherence: classify the nearest pair once per block corner,
    // then shade the block's pixels from its corners' candidate sites by an
    // exact top-2 dot scan. A cell missed by all four corners is dropped.
    const int BLOCK = coherence_block();
    if (BLOCK == 1) {
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
        canvas, BLOCK, positions, scratch_arena_a,
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

  static constexpr int SITES_PER_CORNER = 2;
  /** @brief Current sites endpoint_pairs() finds: the pair and the site that
   *  bounds the others' approach. */
  static constexpr int NEARBY_SITES = SITES_PER_CORNER + 1;
  /** @brief Relative margin on a search bound taken from known sites, for
   *  rounding between the bound's distances and the search's. */
  static constexpr float BOUND_SLACK = 1.0f + 1e-4f;
  /// Corner sites classified at both endpoints, previous first.
  static constexpr int TEMPORAL_SITES_PER_CORNER = 2 * SITES_PER_CORNER;
  /// Margin on the drift bound for rounding in the spin and renormalization.
  static constexpr float DRIFT_SLACK = 1e-5f;

  /** @brief Advances site positions and retains the drawn frame endpoints. */
  HS_HOT_FLASH_MEMBER void advance_sites() {
    HS_PROFILE(vo_animate);
    const bool RESEEDED = active_site_count() != current_num_sites;
    if (RESEEDED)
      seed_sites();

    float s = logf(params.speed + 1.0f) * SITE_SPIN_RADIANS;
    // Every site turns by the same angle, so the half-angle trig is shared;
    // only the axis differs. Expands make_rotation(site.axis, s) exactly.
    const float half_cos = cosf(s / 2);
    const float half_sin = sinf(s / 2);

    previous_positions.clear();
    for (size_t i = 0; i < sites_buffer.size(); ++i) {
      auto &site = sites_buffer[i];
      previous_positions.push_back(site.pos);
      // Renormalize: rotate() drifts |pos| off the unit sphere, and the
      // max-dot nearest search and border acosf(dot) need unit sites.
      const math::Quaternion q(half_cos, half_sin * site.axis);
      site.pos = math::rotate(site.pos, q).normalized();
    }
    frame_spin = s;
    if (RESEEDED || !history_primed) {
      for (size_t i = 0; i < sites_buffer.size(); ++i)
        previous_positions[i] = sites_buffer[i].pos;
      frame_spin = 0.0f;
      history_primed = true;
    }
  }

  /**
   * @brief Coarse-coherence block edge for the active site count.
   * @return Block edge in pixels; it shrinks as cell size falls ~1/sqrt(n).
   */
  int coherence_block() const {
    const float cell_px = (2.0f * H / math::PI_F) /
                          sqrtf(static_cast<float>(sites_buffer.size()));
    return hs::clamp(static_cast<int>(cell_px), COHERENCE_BLOCK_MIN,
                     COHERENCE_BLOCK);
  }

  /** @brief A site's spin from its previous position, as the terms the
   *  spin's cosine and sine scale. */
  struct SiteSpin {
    math::Vector centre; ///< Previous position projected on the spin axis.
    math::Vector cosine; ///< Previous position less its centre.
    math::Vector sine;   ///< Spin axis crossed with the previous position.
  };

  /** @brief Site motion at one column's sweep fraction of the frame. */
  struct ColumnMotion {
    float cosine;          ///< Cosine of the spin so far.
    float sine;            ///< Sine of the spin so far.
    const SiteSpin *spins; ///< Every site's spin terms.
    const KDTree *tree;    ///< Tree over the current positions.
    float drift; ///< Bound on a site's distance from its current position.
    const math::Vector *at; ///< Every site at that time, once filled.
  };

  /**
   * @brief Site motion at a column's sweep time.
   * @param x Absolute canvas column.
   * @param tree Tree over the current positions.
   * @param spins Every site's spin terms.
   * @return The column's motion.
   */
  HS_FLASH_INLINE ColumnMotion column_motion(int x, const KDTree &tree,
                                             const SiteSpin *spins) const {
    const float FRACTION = Temporal::sweep_fraction(x, W);
    const float ANGLE = FRACTION * frame_spin;
    return {cosf(ANGLE),
            sinf(ANGLE),
            spins,
            &tree,
            (1.0f - FRACTION) * fabsf(frame_spin) + DRIFT_SLACK,
            nullptr};
  }

  /**
   * @brief Every site's spin terms for the frame.
   * @param out One entry per site.
   */
  HS_FLASH_INLINE void fill_spins(SiteSpin *out) const {
    for (size_t i = 0; i < sites_buffer.size(); ++i) {
      const math::Vector &from = previous_positions[i];
      const math::Vector &axis = sites_buffer[i].axis;
      const math::Vector centre = axis * math::dot(axis, from);
      out[i] = {centre, from - centre, math::cross(axis, from)};
    }
  }

  /**
   * @brief Fills every site's position at a column's time.
   * @param motion The column's site motion.
   * @param out One position per site.
   */
  __attribute__((always_inline)) void fill_positions(const ColumnMotion &motion,
                                                     math::Vector *out) const {
    for (size_t i = 0; i < sites_buffer.size(); ++i)
      out[i] = position_at(motion, static_cast<uint16_t>(i));
  }

  /**
   * @brief Position of a site part way through the frame's spin.
   * @param motion The column's site motion.
   * @param i Site index.
   * @return The previous position rotated about the site's axis.
   */
  __attribute__((always_inline)) static math::Vector
  position_at(const ColumnMotion &motion, uint16_t i) {
    const SiteSpin &spin = motion.spins[i];
    return spin.centre + spin.cosine * motion.cosine + spin.sine * motion.sine;
  }

  /**
   * @brief Shades a pixel from the exact nearest pair at a column's time.
   * @param p Pixel direction.
   * @param motion The column's site motion; sites move on demand when its
   *        positions are not filled.
   * @param bound_sq Squared distance within which two sites lie at that time,
   *        or FLT_MAX.
   * @return The cell color.
   */
  __attribute__((always_inline)) Color4
  shade_exact(const math::Vector &p, const ColumnMotion &motion,
              float bound_sq = FLT_MAX) const {
    const auto at = [&](uint16_t i) __attribute__((always_inline)) {
      return motion.at ? motion.at[i] : position_at(motion, i);
    };
    const auto nearest = motion.tree->nearest_moved(
        p, 2, motion.drift,
        [&](uint16_t i) { return math::distance_squared(at(i), p); }, bound_sq);
    const uint16_t i0 = nearest[0].original_index;
    const float d0 = math::dot(p, at(i0));
    if (nearest.size() < 2)
      return shade(i0, d0, i0, NO_DOT);
    const uint16_t i1 = nearest[1].original_index;
    return shade(i0, d0, i1, math::dot(p, at(i1)));
  }

  /**
   * @brief Whether a block's corners have the same set of nearest sites at
   *        both endpoints.
   * @param corners The block's corner cells, previous endpoint first.
   * @return False when a site enters or leaves the corners' nearest sites.
   */
  __attribute__((always_inline)) static bool
  same_winners(const std::array<
               const Scan::Shader::BlockCell<TEMPORAL_SITES_PER_CORNER> *, 4>
                   &corners) {
    for (const auto *a : corners) {
      bool before = false, after = false;
      for (const auto *b : corners) {
        before |= (*a)[0] == (*b)[SITES_PER_CORNER];
        after |= (*a)[SITES_PER_CORNER] == (*b)[0];
      }
      if (!before || !after)
        return false;
    }
    return true;
  }

  struct TemporalColumn;

  /** @brief One block's shader at a column's time. */
  struct TemporalBlock {
    TemporalColumn *column;        ///< The column being shaded.
    const math::Vector *positions; ///< Candidate positions at that time.
    const math::Vector *before;    ///< Candidate positions at the frame start.
    const math::Vector *after;     ///< Candidate positions at the frame end.
    const uint16_t *idx;           ///< Candidate site indices.
    uint8_t n;                     ///< Candidate count.
    bool check;                    ///< Whether endpoint winners may differ.

    /** @brief Shades one pixel of the block.
     *  @param p Pixel direction.
     *  @return The cell color. */
    __attribute__((always_inline)) Color4 operator()(const math::Vector &p) {
      float d0 = NO_DOT, d1 = NO_DOT;
      uint8_t b0 = 0, b1 = 0;
      for (uint8_t i = 0; i < n; ++i) {
        const float d = math::dot(p, positions[i]);
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
      if (check && winner(p, before) != winner(p, after))
        return column->effect->shade_exact(
            p, column->motion,
            n > 1 ? fmaxf(math::distance_squared(positions[b0], p),
                          math::distance_squared(positions[b1], p)) *
                        BOUND_SLACK
                  : FLT_MAX);
      return column->effect->shade(idx[b0], d0, idx[b1], d1);
    }

    /** @brief Nearest candidate at one endpoint.
     *  @param p Pixel direction.
     *  @param at Candidate positions at that endpoint.
     *  @return Candidate slot. */
    __attribute__((always_inline)) uint8_t
    winner(const math::Vector &p, const math::Vector *at) const {
      float best = NO_DOT;
      uint8_t slot = 0;
      for (uint8_t i = 0; i < n; ++i) {
        const float d = math::dot(p, at[i]);
        if (d > best) {
          best = d;
          slot = i;
        }
      }
      return slot;
    }
  };

  /** @brief One column's block preparation. */
  struct TemporalColumn {
    const Voronoi *effect;   ///< Sites and shading.
    ColumnMotion motion;     ///< The column's site motion.
    math::Vector *positions; ///< Scratch for one block's candidates, three
                             ///< times its size.

    /**
     * @brief Prepares a block's candidates at the column's time.
     * @param sites The block's corner cells and their union.
     * @return The block's shader, which searches every site where a pixel's
     * nearest candidate differs between the endpoints.
     */
    __attribute__((always_inline)) TemporalBlock
    block(const Scan::Shader::BlockSites<TEMPORAL_SITES_PER_CORNER> &sites) {
      constexpr int SLOTS = 4 * TEMPORAL_SITES_PER_CORNER;
      const bool CHECK = !same_winners(sites.corners);
      math::Vector *before = positions + SLOTS;
      math::Vector *after = before + SLOTS;
      for (uint8_t i = 0; i < sites.n; ++i) {
        const uint16_t SITE = sites.idx[i];
        positions[i] = position_at(motion, SITE);
        if (CHECK) {
          before[i] = effect->previous_positions[SITE];
          after[i] = effect->sites_buffer[SITE].pos;
        }
      }
      return {this, positions, before, after, sites.idx, sites.n, CHECK};
    }
  };

  /**
   * @brief Shades each column at its sweep time between the endpoints.
   * @param canvas Destination canvas.
   * @param tree Tree over the current positions.
   */
  HS_HOT_FLASH_MEMBER void draw_temporal(Canvas &canvas, const KDTree &tree) {
    HS_PROFILE(vo_shade);
    SiteSpin *spins = scratch_arena_a.allocate_n<SiteSpin>(sites_buffer.size());
    fill_spins(spins);
    const int BLOCK = coherence_block();
    if (BLOCK == 1) {
      corners_valid = false;
      math::Vector *every =
          scratch_arena_a.allocate_n<math::Vector>(sites_buffer.size());
      Scan::Shader::draw_columns<W, H>(
          canvas, [&](int x) __attribute__((always_inline)) {
            ColumnMotion motion = column_motion(x, tree, spins);
            fill_positions(motion, every);
            motion.at = every;
            return [this, motion](const math::Vector &p)
                       __attribute__((always_inline)) {
                         return shade_exact(p, motion);
                       };
          });
      return;
    }
    math::Vector *positions = scratch_arena_a.allocate_n<math::Vector>(
        3 * 4 * TEMPORAL_SITES_PER_CORNER);
    Nearby near = first_nearby();
    const bool REUSE =
        corners_valid && corner_clip == canvas.clip() && corner_block == BLOCK;
    size_t corner = 0;
    Scan::Shader::draw_block_columns<W, H, TEMPORAL_SITES_PER_CORNER>(
        canvas, BLOCK, scratch_arena_a,
        [&](const math::Vector &p) __attribute__((always_inline)) {
          return cached_endpoint_pairs(tree, p, near, corner++, REUSE);
        },
        [&](int x) __attribute__((always_inline)) {
          return TemporalColumn{this, column_motion(x, tree, spins), positions};
        });
    corner_clip = canvas.clip();
    corner_block = BLOCK;
    corners_valid = true;
  }

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

  /** @brief Sites nearest a recent sample point, which bound the next
   *  endpoint_pairs() search. */
  using Nearby = std::array<uint16_t, NEARBY_SITES>;

  /**
   * @brief Distinct sites to bound the first endpoint_pairs() search.
   * @return The first NEARBY_SITES sites, repeating the last when fewer exist.
   */
  Nearby first_nearby() const {
    Nearby near{};
    for (int i = 0; i < NEARBY_SITES; ++i)
      near[i] = static_cast<uint16_t>(
          std::min(i, static_cast<int>(sites_buffer.size()) - 1));
    return near;
  }

  /**
   * @brief Nearest pairs of a sample point at both endpoints of the frame's
   *        spin.
   * @param current Tree over the current positions.
   * @param p Sample point on the unit sphere.
   * @param near Sites nearest a recent sample point, which bound the search;
   *        receives this point's nearest current sites.
   * @return The previous endpoint's nearest pair, then the current's, each
   *         nearest first, as an exhaustive search ranks them.
   * @details The current pair is the previous pair too, reordered, when no
   * other site can come within its previous distances after moving the
   * frame's spin; otherwise the tree is searched for the previous positions.
   */
  __attribute__((always_inline))
  Scan::Shader::BlockCell<TEMPORAL_SITES_PER_CORNER>
  endpoint_pairs(const KDTree &current, const math::Vector &p,
                 Nearby &near) const {
    float bound = 0.0f;
    for (uint16_t site : near)
      bound = fmaxf(bound, math::distance_squared(sites_buffer[site].pos, p));
    const auto knn = current.nearest_moved(
        p, NEARBY_SITES, 0.0f,
        [&](uint16_t i) {
          return math::distance_squared(sites_buffer[i].pos, p);
        },
        bound * BOUND_SLACK);
    for (size_t i = 0; i < near.size(); ++i)
      near[i] = knn[std::min(i, knn.size() - 1)].original_index;
    const uint16_t a = near[0];
    const uint16_t b = near[1];
    const float from_a = math::distance_squared(previous_positions[a], p);
    const float from_b = math::distance_squared(previous_positions[b], p);
    if (knn.size() > 2) {
      const float spin = fabsf(frame_spin) + DRIFT_SLACK;
      const float reach = sqrtf(knn[2].d_sq) - spin;
      const float from = fmaxf(from_a, from_b) * BOUND_SLACK;
      if (!(reach > 0.0f && from < reach * reach)) {
        const auto pair = current.nearest_moved(
            p, SITES_PER_CORNER, spin,
            [&](uint16_t i) {
              return math::distance_squared(previous_positions[i], p);
            },
            from);
        return {pair[0].original_index, pair[1].original_index, a, b};
      }
    }
    const bool b_first = from_b < from_a || (from_b == from_a && b < a);
    return {b_first ? b : a, b_first ? a : b, a, b};
  }

  /** @brief Endpoint pair retained in corner traversal order. */
  using CornerPair = std::array<uint16_t, SITES_PER_CORNER>;

  /**
   * @brief Classifies a corner with its retained previous endpoint pair.
   * @param current Tree over current site positions.
   * @param p Corner direction.
   * @param near Search bounds for an uncached traversal.
   * @param corner Corner ordinal in the unchanged clip and grid.
   * @param reuse Whether the preceding frame classified the same grid.
   * @return Previous and current nearest pairs, nearest first.
   */
  HS_HOT_INLINE Scan::Shader::BlockCell<TEMPORAL_SITES_PER_CORNER>
  cached_endpoint_pairs(const KDTree &current, const math::Vector &p,
                        Nearby &near, size_t corner, bool reuse) {
    auto &saved = corner_pairs[corner];
    if (!reuse) {
      const auto pairs = endpoint_pairs(current, p, near);
      saved = {pairs[2], pairs[3]};
      return pairs;
    }
    const CornerPair before = saved;
    const float BOUND =
        fmaxf(math::distance_squared(sites_buffer[before[0]].pos, p),
              math::distance_squared(sites_buffer[before[1]].pos, p));
    const auto pair = current.nearest_moved(
        p, SITES_PER_CORNER, 0.0f,
        [&](uint16_t i) HS_HOT_INLINE {
          return math::distance_squared(sites_buffer[i].pos, p);
        },
        BOUND * BOUND_SLACK);
    saved = {pair[0].original_index,
             pair[pair.size() > 1 ? 1 : 0].original_index};
    return {before[0], before[1], saved[0], saved[1]};
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

  using CandSet = Scan::Shader::BlockCandidates<SITES_PER_CORNER>;

  // Compile-time high-water check for the 64 KB scratch_arena_a reserve. The
  // transient peaks share that arena, with positions + KD nodes live across all:
  //   build:    positions + KD nodes + KD build-index scratch
  //   shading:  positions + KD nodes + corner cells + candidate-row sets
  //   temporal: positions + KD nodes + corner cells + block-column
  //             site unions + candidate and per-site positions + spin terms
  static constexpr size_t SCRATCH_A_BYTES = 64 * 1024;
  static constexpr size_t POSITIONS_BYTES =
      size_t(MAX_SITES) * sizeof(math::Vector);
  static constexpr size_t KD_NODES_BYTES = size_t(MAX_SITES) * sizeof(KDNode);
  static constexpr size_t KD_BUILD_SCRATCH_BYTES =
      size_t(MAX_SITES) * sizeof(int);
  static constexpr size_t WALKER_BYTES =
      Scan::Shader::block_coherent_scratch_bytes<W, H, SITES_PER_CORNER,
                                                 COHERENCE_BLOCK_MIN>();
  static constexpr size_t TEMPORAL_WALKER_BYTES =
      Scan::Shader::block_columns_scratch_bytes<
          W, H, TEMPORAL_SITES_PER_CORNER, std::max(2, COHERENCE_BLOCK_MIN)>() +
      (3 * 4 * TEMPORAL_SITES_PER_CORNER + MAX_SITES) * sizeof(math::Vector) +
      MAX_SITES * sizeof(SiteSpin);
  static constexpr size_t SCRATCH_HIGH_WATER =
      std::max(POSITIONS_BYTES + KD_NODES_BYTES +
                   std::max(KD_BUILD_SCRATCH_BYTES, WALKER_BYTES),
               POSITIONS_BYTES + KD_NODES_BYTES +
                   std::max(KD_BUILD_SCRATCH_BYTES, TEMPORAL_WALKER_BYTES));
  static_assert(
      SCRATCH_HIGH_WATER <= SCRATCH_A_BYTES,
      "Voronoi scratch_arena_a budget too small for MAX_SITES positions + "
      "KD-tree + coarse-grid cells; raise SCRATCH_A_BYTES or lower MAX_SITES / "
      "coarsen COHERENCE_BLOCK");

  /// Maximum corners in a block grid that uses temporal coherence.
  static constexpr size_t CORNER_CAPACITY =
      size_t((W - 1) / std::max(2, COHERENCE_BLOCK_MIN) + 2) *
      size_t((H - 1) / std::max(2, COHERENCE_BLOCK_MIN) + 2);

  // The sites and their previous positions live in the persistent arena.
  static constexpr size_t FOOTPRINT_BYTES =
      size_t(MAX_SITES) * sizeof(Site) + alignof(Site) +
      size_t(MAX_SITES) * sizeof(math::Vector) + alignof(math::Vector) +
      CORNER_CAPACITY * sizeof(CornerPair) + alignof(CornerPair);
  static_assert(FOOTPRINT_BYTES <=
                    ArenaSplit{SCRATCH_A_BYTES, 0}.device_persistent(),
                "Voronoi persistent footprint exceeds its device partition; "
                "lower MAX_SITES or shrink SCRATCH_A_BYTES");

  int current_num_sites = 0;      /**< Count currently seeded; re-seeds (clear +
                                   refill, no realloc) when the slider changes. */
  ArenaVector<Site> sites_buffer; /**< Active Voronoi sites for the frame. */
  /// Site positions at the previous frame's end, parallel to `sites_buffer`.
  ArenaVector<math::Vector> previous_positions;
  float frame_spin = 0.0f;     ///< Spin from the previous positions, radians.
  bool history_primed = false; ///< Whether a previous endpoint was drawn.

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

  CornerPair *corner_pairs = nullptr; ///< Current endpoint pairs by corner.
  ClipRegion corner_clip;     ///< Clip that produced the retained pairs.
  int corner_block = 0;       ///< Retained grid's block edge.
  bool corners_valid = false; ///< Pairs belong to the preceding draw.
};
