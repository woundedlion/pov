/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file spherical_tiling.h
 * @brief Shader interpreter crossing for registry mesh tiles. */

#include "core/platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include <array>
#include <memory>

#include "core/math/scalar_kernels.h"
#include "core/pattern/field.h"
#include "core/render/pullback/operators/common.h"
#include "core/render/pullback/sphere_tiling.h"
#include "core/render/pullback/spherical_tiling_source.h"

namespace Pullback::Interp::Op {

struct SphericalTilingGeometryCache {
  std::array<std::shared_ptr<const Source::SphereTilingGeometry>, 2> geometry{};
  std::array<uint8_t, 2> solids{};
};

/** @brief Projected geometry for both endpoints of a source morph. */
struct SphericalTilingState {
  mutable std::shared_ptr<const SphericalTilingGeometryCache> cache;
  mutable uint8_t active = 0;
};

/** @brief Sphere-to-field crossing for tile fills and seams. */
struct SampleSphericalTiling : ValueStateModel<SphericalTilingState> {
  static constexpr const char *ID =
      Source::SphericalTilingSourceParams::OPERATOR_ID;
  static constexpr const char *NAME = "Tiling";
  using Input = SphereSample;
  using Output = FieldSample;
  using Params = Source::SphericalTilingSourceParams;

  /** @brief Frame controls and a view owned by the operator state. */
  struct Prepared {
    const Source::SphereTilingGeometry *geometry; ///< Selected solid index.
    Source::SphereTilingEdgeBand edge;            ///< Edge coverage controls.
    float half_width;                             ///< Face inset in radians.
    float fringe;                                 ///< Face depth fraction.
  };

  static void advance(State &, const Params &) {}

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    const auto &geometry = geometry_for(state, params.solid);
    const Source::SphereTilingEdgeBand edge{params.line_thickness,
                                            params.line_softness};
    return {&geometry, edge, edge.half_width,
            params.line_softness / Params::MAX_SOFTNESS};
  }

  static FieldSample run(const SphereSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const auto direction = Pattern::SphereDirection::admit(input.dir);
    const auto &geometry = *prepared.geometry;
    const auto output = static_cast<Source::TilingOutput>(params.output);
    if (output == Source::TilingOutput::EDGES) {
      const auto edge =
          geometry.edge_view().sample(direction.value(), prepared.edge);
      return {edge.value, edge.coverage, input.dir, input.path_length};
    }
    const uint16_t face = geometry.face_at(direction.value());
    if (face == Source::SphereTilingGeometry::NO_FACE)
      return {0.0f, 0.0f, input.dir, input.path_length};
    const float radius = ctx.angular_sample_radius;
    const auto depth = geometry.face_depth(face, direction.value(), radius);
    const float span = depth.center - prepared.half_width;
    const float position = hs::clamp((depth.sample - prepared.half_width) /
                                         std::max(math::TOLERANCE, span),
                                     0.0f, 1.0f);
    const float softness = std::max(0.0f, span) * prepared.fringe;
    const float across =
        std::max(0.0f, depth.across_seam - prepared.half_width) *
        prepared.fringe;
    const float alpha = face_coverage(depth.sample, radius, prepared.half_width,
                                      softness, across);
    return {position, alpha, input.dir, input.path_length};
  }

private:
  static const Source::SphereTilingGeometry &geometry_for(const State &state,
                                                          uint8_t solid) {
    if (state.cache)
      for (uint8_t index = 0; index < state.cache->geometry.size(); ++index)
        if (state.cache->geometry[index] &&
            state.cache->solids[index] == solid) {
          state.active = index;
          return *state.cache->geometry[index];
        }
    auto cache =
        state.cache
            ? std::make_shared<SphericalTilingGeometryCache>(*state.cache)
            : std::make_shared<SphericalTilingGeometryCache>();
    state.active ^= 1;
    cache->geometry[state.active] =
        std::make_shared<const Source::SphereTilingGeometry>(solid);
    cache->solids[state.active] = solid;
    state.cache = std::move(cache);
    return *state.cache->geometry[state.active];
  }

  static float integrated_fade(float depth, float softness) {
    if (depth <= 0.0f)
      return 0.0f;
    if (depth >= softness)
      return depth - 0.5f * softness;
    return 0.5f * depth * (depth / softness);
  }

  static float face_coverage(float distance, float radius, float half_width,
                             float softness, float across_softness) {
    const float inner = distance - radius;
    const float outer = distance + radius;
    if (outer <= half_width)
      return 0.0f;
    if (inner - half_width >= softness)
      return 1.0f;
    if (radius == 0.0f)
      return (distance - half_width) / softness;
    const float outer_area = integrated_fade(outer - half_width, softness);
    const float inner_area = integrated_fade(
        fabsf(inner) - half_width, inner < 0.0f ? across_softness : softness);
    // The footprint can include the inward fade on both sides of the seam.
    const float area = outer_area - copysignf(inner_area, inner);
    return hs::clamp(area / (2.0f * radius), 0.0f, 1.0f);
  }
};

} // namespace Pullback::Interp::Op

#endif // HS_ENABLE_CHAIN_INTERPRETER
