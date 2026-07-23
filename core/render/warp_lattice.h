/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "engine/memory.h"
#include "engine/platform.h"
#include "engine/reaction_graph.h"
#include "math/3dmath.h"

/**
 * @brief Sphere-uniform vector field sampled on the ReactionGraph Fibonacci
 *        lattice: populate() evaluates a warp at every node, sample()
 *        reconstructs it at any direction.
 * @details An equirect coarse grid cannot represent a smooth sphere warp near
 * the poles — in map coordinates the field's gradient grows as 1/sin(phi) and
 * pole-crossing targets jump by half the map width. The lattice samples and
 * interpolates entirely in 3D, where the field is smooth everywhere, so
 * reconstruction quality is latitude-independent. Node positions are stored as
 * snorm16 x/z pairs (y is analytic in the node index) and targets as snorm16
 * triples; kernel weights mirror ReactionDiffusionBase's Wendland C2 stencil
 * (effects/ReactionDiffusionBase.h). The final renormalization absorbs the
 * kernel's weight sum, so sample() needs no divide.
 */
class WarpLattice {
public:
  static constexpr int N = ReactionGraph::RD_N;
  static constexpr int K = ReactionGraph::RD_K;

  /** @brief Kernel support radius over the mean node spacing. */
  static constexpr float KERNEL_R = 1.5f * ReactionGraph::D_AVG;
  static constexpr float INV_R2 = 1.0f / (KERNEL_R * KERNEL_R);

  /** @brief Contiguous node-index range; the lattice is ordered north (i = 0)
   *  to south, so a latitude band is an index interval. */
  struct Band {
    int lo, hi;
    bool operator==(const Band &) const = default;
  };

  /** @brief Persistent bytes init_storage() reserves (cubemap LUT + node
   *  positions + targets). The cubemap build additionally carves a ~90 KB
   *  transient that is rewound before the arrays are placed. */
  static constexpr size_t STORAGE_BYTES =
      6u * ReactionGraph::CubemapLUT::RES * ReactionGraph::CubemapLUT::RES *
          sizeof(uint16_t) +
      2u * N * sizeof(int16_t) + 3u * N * sizeof(int16_t);

  /**
   * @brief Allocates the cubemap LUT and node/target arrays from @p arena.
   * @param arena Persistent arena; must have STORAGE_BYTES plus the cubemap
   *        build's transient headroom available.
   * @details Call from effect init() (arenas are not ready in constructors),
   * and again after any compaction that resets the arena.
   */
  HS_COLD_MEMBER void init_storage(Arena &arena) {
    cube_lut.build(arena);
    pos_xz = arena.allocate_n<int16_t>(2 * N);
    targets = arena.allocate_n<int16_t>(3 * N);
    for (int i = 0; i < N; ++i) {
      Vector p = ReactionGraph::node(i);
      pos_xz[2 * i] = quantize(p.x);
      pos_xz[2 * i + 1] = quantize(p.z);
    }
  }

  /** @brief True once init_storage() has bound the arrays. */
  bool bound() const { return targets != nullptr; }

  /** @brief 3D y of lattice node @p i (exact, no table). */
  HS_O3_FN static float node_y(int i) {
    return 1.0f - 2.0f * static_cast<float>(i) / (N - 1);
  }

  /**
   * @brief Node-index band covering the 3D y range [y_south, y_north],
   *        widened by the kernel's gather reach.
   * @param y_north Larger 3D y (closer to the north pole, smaller index).
   * @param y_south Smaller 3D y.
   */
  static Band band(float y_north, float y_south) {
    // Index margin for one kernel radius plus a neighbor ring: |dy| <= angle,
    // and index advances (N - 1) / 2 per unit y.
    constexpr int MARGIN =
        static_cast<int>(2.0f * ReactionGraph::D_AVG * 0.5f * (N - 1)) + 1;
    int lo = static_cast<int>((1.0f - y_north) * 0.5f * (N - 1)) - MARGIN;
    int hi = static_cast<int>((1.0f - y_south) * 0.5f * (N - 1)) + MARGIN + 1;
    return {hs::clamp(lo, 0, N), hs::clamp(hi, 0, N)};
  }

  /**
   * @brief Evaluates @p warp at every node in @p b and stores the targets.
   * @tparam WarpFn Callable mapping a unit Vector to a unit Vector.
   * @param b Node range to populate (see band()).
   * @param warp The field; evaluated at the quantized node position so
   *        populate and sample agree on geometry exactly.
   */
  template <typename WarpFn> void populate(Band b, WarpFn &&warp) {
    for (int i = b.lo; i < b.hi; ++i) {
      Vector t = warp(position(i));
      targets[3 * i] = quantize(t.x);
      targets[3 * i + 1] = quantize(t.y);
      targets[3 * i + 2] = quantize(t.z);
    }
  }

  /**
   * @brief Reconstructs the field at direction @p v.
   * @param v Unit query direction.
   * @return Unit target direction: the Wendland-weighted blend of the nearby
   *         nodes' targets, renormalized.
   * @details The cubemap seed can be off by a neighbor, so the stencil tracks
   * its argmin and re-walks from the genuine nearest node when the seed loses —
   * the same bias removal as ReactionDiffusionBase::refine_and_accumulate. A
   * near-zero blend (cancelling targets) falls back to the nearest node's
   * target rather than normalizing a degenerate vector.
   */
  HS_O3_FN Vector sample(const Vector &v) const {
    float d2s[K + 1];
    int ids[K + 1];
    int n = stencil(cube_lut.lookup(v), v, d2s, ids);
    int best = 0;
    for (int i = 1; i < n; ++i)
      if (d2s[i] < d2s[best]) best = i;
    if (best != 0) {
      n = stencil(ids[best], v, d2s, ids);
      best = 0;
    }
    Vector acc(0.0f, 0.0f, 0.0f);
    for (int i = 0; i < n; ++i) {
      float u = 1.0f - d2s[i] * INV_R2;
      if (u > 0.0f) acc += target(ids[i]) * (u * u);
    }
    if (dot(acc, acc) < 1e-12f) return target(ids[best]);
    return acc.normalized();
  }

private:
  static int16_t quantize(float c) {
    return static_cast<int16_t>(hs::clamp(c, -1.0f, 1.0f) * 32767.0f);
  }
  static constexpr float DEQUANT = 1.0f / 32767.0f;

  HS_O3_FN Vector position(int i) const {
    return Vector(pos_xz[2 * i] * DEQUANT, node_y(i),
                  pos_xz[2 * i + 1] * DEQUANT);
  }

  HS_O3_FN Vector target(int i) const {
    return Vector(targets[3 * i] * DEQUANT, targets[3 * i + 1] * DEQUANT,
                  targets[3 * i + 2] * DEQUANT);
  }

  /** @brief Fills @p d2s / @p ids with @p center and its neighbors' squared
   *  distances to @p v; returns the entry count. */
  HS_O3_FN int stencil(int center, const Vector &v, float *d2s,
                       int *ids) const {
    auto d2 = [&](int node) {
      Vector p = position(node);
      float dx = v.x - p.x, dy = v.y - p.y, dz = v.z - p.z;
      return dx * dx + dy * dy + dz * dz;
    };
    d2s[0] = d2(center);
    ids[0] = center;
    int n = 1;
    for (int k = 0; k < K; ++k) {
      int ni = ReactionGraph::neighbors[center][k];
      if (ni < 0) continue;
      d2s[n] = d2(ni);
      ids[n] = ni;
      ++n;
    }
    return n;
  }

  ReactionGraph::CubemapLUT cube_lut;
  int16_t *pos_xz = nullptr;
  int16_t *targets = nullptr;
};
