/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file spherical_tiling.h
 * @brief Shader interpreter crossing for registry mesh tiles. */

#include "core/platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include <memory>

#include "core/math/scalar_kernels.h"
#include "core/pattern/field.h"
#include "core/render/pullback/operators/common.h"
#include "core/render/pullback/sphere_tiling.h"
#include "core/render/pullback/spherical_tiling_source.h"

namespace Pullback::Interp::Op {

/** @brief Per-instance cache of projected geometry for the selected solid. */
struct SphericalTilingState {
  mutable std::shared_ptr<const Source::SphereTilingGeometry>
      geometry;              ///< Immutable index shared by migrated instances.
  mutable uint8_t solid = 0; ///< Registry ordinal stored in geometry.
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
    float half_width; ///< Edge half-width in radians.
    float fringe_end; ///< Outer limit of the opacity fringe in radians.
  };

  static void advance(State &, const Params &) {}

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    if (!state.geometry || state.solid != params.solid) {
      state.geometry =
          std::make_shared<const Source::SphereTilingGeometry>(params.solid);
      state.solid = params.solid;
    }
    constexpr float DEGREES_TO_RADIANS = math::PI_F / 180.0f;
    const float half_width = 0.5f * params.line_thickness;
    return {state.geometry.get(), half_width * DEGREES_TO_RADIANS,
            (half_width + params.line_softness) * DEGREES_TO_RADIANS};
  }

  static FieldSample run(const SphereSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    const auto direction = Pattern::SphereDirection::admit(input.dir);
    const auto &geometry = *prepared.geometry;
    const auto output = static_cast<Source::TilingOutput>(params.output);
    if (output == Source::TilingOutput::EDGES) {
      const auto edge =
          geometry.edge_at(direction.value(), prepared.fringe_end);
      const float alpha =
          1.0f - math::smooth_ramp(prepared.half_width, prepared.fringe_end,
                                   edge.distance);
      return {edge.position, alpha, input.dir, input.path_length};
    }
    const uint16_t face = geometry.face_at(direction.value());
    if (face == Source::SphereTilingGeometry::NO_FACE)
      return {0.0f, 0.0f, input.dir, input.path_length};
    const auto depth = geometry.face_depth(face, direction.value());
    const float span =
        std::max(math::TOLERANCE, depth.center - prepared.half_width);
    const float position =
        hs::clamp((depth.sample - prepared.half_width) / span, 0.0f, 1.0f);
    const float alpha = math::smooth_ramp(prepared.half_width,
                                          prepared.fringe_end, depth.sample);
    return {position, alpha, input.dir, input.path_length};
  }
};

} // namespace Pullback::Interp::Op

#endif // HS_ENABLE_CHAIN_INTERPRETER
