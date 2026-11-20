/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/pullback/composed_effect.h"
#include "core/render/pullback/sphere_tiling_edges_baked.h"
#include "core/render/pullback/spherical_tiling_source.h"

/** @file composed_tiling.h
 * @brief Compiled spherical tiling edge source.
 */

namespace Pullback {

namespace Source {

inline bool valid(const SphericalTilingSourceParams &params) {
  return Fields::valid(params) &&
         params.output == static_cast<uint8_t>(TilingOutput::EDGES) &&
         supports_baked_sphere_tiling_edges(params.solid);
}

/** @brief Samples the nearest edge in a baked spherical tiling.
 * @tparam State Provider of tiling parameters.
 */
template <typename State> struct BakedTilingEdges : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  struct Prepared {
    SphereTilingEdgeView edges;
    SphereTilingEdgeBand band{0.0f, 0.05f};
  };

  static Prepared prepare(const FrameState &frame) {
    const auto &p = State::params(frame);
    return {baked_sphere_tiling_edges(p.solid),
            SphereTilingEdgeBand(p.line_thickness, p.line_softness)};
  }

  static HS_HOT_FLASH_MEMBER FieldSample sample(const SphereSample &input,
                                                const FrameState &,
                                                const Prepared &source) {
    const auto direction = Pattern::SphereDirection::admit(input.dir);
    const auto edge = source.edges.sample(direction.value(), source.band);
    return {edge.value, edge.coverage, input.dir, input.path_length};
  }
};

} // namespace Source

/** @brief Baked tiling edge source. @tparam B Frame binding. */
template <typename B>
struct SourcePolicyFor<Source::SphericalTilingSourceParams, B> {
  using Type = Source::BakedTilingEdges<
      SourceProvider<B, Source::SphericalTilingSourceParams>>;
};

} // namespace Pullback
