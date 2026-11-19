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
#include "core/math/scalar_noise.h"
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
    float border_start; ///< Border half-width in radians.
    float border_end;   ///< Outer fade limit in radians.
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
    return {state.geometry.get(), params.line_thickness * DEGREES_TO_RADIANS,
            (params.line_thickness + params.line_softness) *
                DEGREES_TO_RADIANS};
  }

  static FieldSample run(const SphereSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    const auto direction = Pattern::SphereDirection::admit(input.dir);
    const auto &geometry = *prepared.geometry;
    const auto output = static_cast<Source::TilingOutput>(params.output);
    if (output == Source::TilingOutput::BORDERS) {
      const float distance =
          geometry.border_distance(direction.value(), prepared.border_end);
      const float shade = math::smooth_ramp(prepared.border_start,
                                            prepared.border_end, distance);
      return Kernel::sample(input, 1.0f - 2.0f * shade);
    }
    const uint16_t face = geometry.face_at(direction.value());
    if (face == Source::SphereTilingGeometry::NO_FACE)
      return Kernel::sample(input, -1.0f);
    const uint16_t value = output == Source::TilingOutput::FACE_TYPE
                               ? geometry.face_type(face)
                               : face;
    return Kernel::sample(input, 2.0f * math::hash01(value, 0u) - 1.0f);
  }
};

} // namespace Pullback::Interp::Op

#endif // HS_ENABLE_CHAIN_INTERPRETER
