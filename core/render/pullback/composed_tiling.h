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

  struct Endpoint {
    SphereTilingEdgeView edges;
    SphereTilingEdgeBand band{0.0f, 0.05f};
  };

  struct Morph {
    Endpoint incoming;
    float mix = 0.0f;
  };

  struct NoMorph {};

  static constexpr bool MORPHS =
      requires(const FrameState &frame) { frame.source_morph; };

  struct Prepared {
    Endpoint first;
#if defined(_MSC_VER)
    [[msvc::no_unique_address]]
#else
    [[no_unique_address]]
#endif
    std::conditional_t<MORPHS, Morph, NoMorph> morph{};
  };

  static Prepared prepare(const FrameState &frame) {
    if constexpr (MORPHS) {
      const auto &morph = frame.source_morph;
      if (morph.from && morph.fields_differ) {
        if (morph.mix <= 0.0f)
          return {prepare_endpoint(morph.from->template get<"source">())};
        auto incoming = prepare_endpoint(morph.to->template get<"source">());
        if (morph.mix >= 1.0f)
          return {incoming};
        return {prepare_endpoint(morph.from->template get<"source">()),
                {incoming, morph.mix}};
      }
    }
    return {prepare_endpoint(State::params(frame))};
  }

  static HS_HOT_INLINE FieldSample sample(const SphereSample &input,
                                          const FrameState &,
                                          const Prepared &source) {
    const auto direction = Pattern::SphereDirection::admit(input.dir);
    const auto projection = math::project_cube_face(direction.value());
    const int narrow_cell = SphereTilingNarrowIndex::cell_for(projection);
    auto edge = sample_endpoint(direction.value(), source.first, projection,
                                narrow_cell);
    if constexpr (MORPHS) {
      if (source.morph.mix > 0.0f) {
        const auto incoming = sample_endpoint(
            direction.value(), source.morph.incoming, projection, narrow_cell);
        edge = blend_tiling_edges(edge, incoming, source.morph.mix);
      }
    }
    return {edge.value, edge.coverage, input.dir, input.path_length};
  }

private:
  static Endpoint prepare_endpoint(const SphericalTilingSourceParams &params) {
    return {baked_sphere_tiling_edges(params.solid),
            SphereTilingEdgeBand(params.line_thickness, params.line_softness)};
  }

  static HS_HOT_ITCM_MEMBER SphereTilingEdgeValue
  sample_endpoint(const math::Vector &direction, const Endpoint &endpoint,
                  const math::CubeFaceProjection &projection, int narrow_cell) {
    return endpoint.edges.sample(direction, endpoint.band, projection,
                                 narrow_cell);
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
