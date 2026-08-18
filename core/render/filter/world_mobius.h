/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/filter/pipeline.h"
#include "math/stereographic.h"
#include "math/geometry.h"
#include "color/color.h"

/**
 * @file world_mobius.h
 * @brief Filter::World::Mobius: stereographically projects, applies a live
 * Mobius map, and re-emits on the sphere.
 */

namespace Filter {
namespace World {

/**
 * @brief Applies a Mobius transformation to 3D points.
 */
class Mobius : public Is3D {
public:
  static constexpr bool requires_unit_world_input = true;
  /**
   * @brief The map is non-rigid, so no rotation-mirroring cull_edge can bound
   *        the image of an edge; the effect must render the full canvas.
   */
  static constexpr bool crosses_segments = true;
  /**
   * @brief Binds the filter to a live Mobius parameter set.
   * @param params Mobius transform parameters applied per point.
   */
  Mobius(MobiusParams &params) : params(params) {}
  /**
   * @brief Stereographically projects, applies the Mobius map, and re-emits.
   * @param v World-space point on the unit sphere.
   * @param color Source color, forwarded unchanged.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    pass(inv_stereo(mobius(stereo(v), params)), color, age, alpha);
  }

private:
  MobiusParams &params; /**< Live Mobius transform parameters. */
};

} // namespace World
} // namespace Filter
