/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/filter/pipeline.h"
#include "math/geometry.h"
#include "color/color.h"

/**
 * @file world_replicate.h
 * @brief Filter::World::Replicate: emits evenly-spaced copies of the geometry
 * rotated around the Y axis.
 */

namespace Filter {
namespace World {

/**
 * @brief Replicates geometry by rotating it around the Y-axis.
 * @details Every copy shares the source `age` (replication is spatial, not
 * temporal).
 */
template <int W> class Replicate : public Is3D {
public:
  /**
   * @brief Builds a replicator emitting @p count evenly-spaced Y-axis copies.
   * @param count Desired copy count; clamped to [1, W].
   * @details `this->count` (the clamped member, declared/initialized first)
   * feeds make_rotation, so count == 0 cannot feed inf into it. The ceiling is
   * W because W copies already sit one equatorial pixel column apart; beyond
   * that they land on the same column.
   */
  Replicate(int count) { set_count(count); }

  /**
   * @brief Changes the number of evenly spaced copies.
   * @param new_count Desired copy count; clamped to [1, W].
   */
  void set_count(int new_count) {
    count = hs::clamp(new_count, 1, W);
    step = make_rotation(Y_AXIS, 2 * PI_F / count);
  }
  /**
   * @brief Emits the point plus count-1 rotated copies around the Y axis.
   * @param v World-space point to replicate.
   * @param color Source color, forwarded unchanged to every copy.
   * @param age Temporal age channel (frames), shared by every copy.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    Vector r = v;
    pass(r, color, age, alpha);
    for (int i = 1; i < count; i++) {
      // renormalize so repeated rotation can't drift copies off the unit sphere
      r = rotate(r, step).normalized();
      pass(r, color, age, alpha);
    }
  }

  /**
   * @brief Re-emits a clip-cull edge under every copy's Y-axis rotation.
   * @tparam FwdFn Downstream cull continuation
   *         `bool(const Vector&, const Vector&, const Basis*)`.
   * @param a,b Edge endpoints in world space (pre-rotation).
   * @param pb Optional planar basis, rotated alongside the endpoints.
   * @param forward Tail-of-pipeline cull continuation.
   * @return True if any copy of the edge could intersect the clip band.
   * @details Mirrors plot()'s rotation loop so the cull sees the longitudes the
   *          copies are drawn at, not the source geometry's.
   */
  template <typename FwdFn>
  bool cull_edge(const Vector &a, const Vector &b, const Basis *pb,
                 FwdFn &&forward) const {
    if (forward(a, b, pb))
      return true;
    Vector ra = a, rb = b;
    Basis rp;
    if (pb)
      rp = *pb;
    for (int i = 1; i < count; i++) {
      ra = rotate(ra, step).normalized();
      rb = rotate(rb, step).normalized();
      if (pb) {
        rp = rotate(rp, step);
        if (forward(ra, rb, &rp))
          return true;
      } else if (forward(ra, rb, nullptr)) {
        return true;
      }
    }
    return false;
  }

private:
  int count;       /**< Number of copies emitted, in [1, W]. */
  Quaternion step; /**< Per-copy Y-axis rotation (2*pi / count). */
};

} // namespace World
} // namespace Filter
