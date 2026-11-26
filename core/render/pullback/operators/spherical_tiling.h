/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file spherical_tiling.h
 * @brief Shader interpreter crossing for registry mesh tiles. */

#include "core/platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include <algorithm>
#include <array>
#include <memory>

#include "core/math/scalar_kernels.h"
#include "core/pattern/field.h"
#include "core/render/pullback/operators/common.h"
#include "core/render/pullback/sphere_tiling.h"
#include "core/render/pullback/spherical_tiling_source.h"

namespace Pullback::Interp::Op {

/** @brief Heap-owned faces-and-edges geometry of one registry solid. */
class OwnedSphereTiling {
public:
  /** @brief Builds a solid's geometry and keeps it in an exact-size block.
   * @param solid Registry solid ordinal. */
  explicit OwnedSphereTiling(uint8_t solid) {
    constexpr size_t BUILD_BYTES = 3 * 1024 * 1024;
    constexpr size_t SCRATCH_BYTES = 512 * 1024;
    auto build = std::make_unique_for_overwrite<uint8_t[]>(BUILD_BYTES +
                                                           2 * SCRATCH_BYTES);
    const auto scratch_a =
        std::make_unique<Arena>(build.get() + BUILD_BYTES, SCRATCH_BYTES);
    const auto scratch_b = std::make_unique<Arena>(
        build.get() + BUILD_BYTES + SCRATCH_BYTES, SCRATCH_BYTES);
    using Geometry = Source::SphereTilingGeometry;
    const auto built = Geometry::build<Geometry::Content::FACES_AND_EDGES>(
        solid, {build.get(), BUILD_BYTES}, *scratch_a, *scratch_b);
    block = std::make_unique_for_overwrite<uint8_t[]>(built.bytes());
    std::copy_n(build.get(), built.bytes(), block.get());
    geometry = built.rebased(build.get(), block.get());
  }

  /** @brief The solid's geometry.
   * @return Geometry viewing the owned block. */
  const Source::SphereTilingGeometry &get() const { return geometry; }

private:
  std::unique_ptr<uint8_t[]> block;
  Source::SphereTilingGeometry geometry;
};

/// Shared immutable geometry for the two most recently selected distinct solids.
struct SphericalTilingGeometryCache {
  std::array<std::shared_ptr<const OwnedSphereTiling>, 2>
      geometry{}; ///< Projected geometry; null until populated.
  std::array<uint8_t, 2>
      solids{}; ///< Registry solid ordinal for each populated slot.
};

/** @brief Projected geometry for both endpoints of a source morph. */
struct SphericalTilingState {
  mutable std::shared_ptr<const SphericalTilingGeometryCache>
      cache; ///< Shared geometry snapshot; replaced on a cache miss.
  mutable uint8_t active = 0; ///< Most recently selected cache slot.
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
    const Source::SphereTilingGeometry *geometry; ///< Selected solid geometry.
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
          return state.cache->geometry[index]->get();
        }
    auto cache =
        state.cache
            ? std::make_shared<SphericalTilingGeometryCache>(*state.cache)
            : std::make_shared<SphericalTilingGeometryCache>();
    state.active ^= 1;
    cache->geometry[state.active] =
        std::make_shared<const OwnedSphereTiling>(solid);
    cache->solids[state.active] = solid;
    state.cache = std::move(cache);
    return state.cache->geometry[state.active]->get();
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
