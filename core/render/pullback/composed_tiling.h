/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include "core/render/pullback/composed_effect.h"
#include "core/render/pullback/sphere_tiling.h"
#include "core/render/pullback/spherical_tiling_source.h"

/** @file composed_tiling.h
 * @brief Compiled spherical tiling edge source.
 */

namespace Pullback {

namespace Source {

/** @brief Block bytes of each registry solid's larger edge-only geometry.
 * @details Host results of SphereTilingGeometry::build with Content::EDGES,
 * the larger of the Indices::NARROW and Indices::COARSE builds, in flattened
 * registry order; a test rebuilds every entry. */
inline constexpr std::array<uint32_t, Solids::NUM_ENTRIES>
    SPHERE_TILING_EDGE_BYTES = {
        11260,  11956,  11404,  15436,  13516,  14884,  13444,  17668,
        17332,  16276,  21700,  16804,  19180,  26092,  22540,  24964,
        32692,  27364,  12388,  14236,  14476,  14380,  16492,  18508,
        18796,  17428,  19276,  20044,  24028,  28204,  26788,  94516,
        148488, 153504, 206644, 105648, 268900, 107736, 223780, 203044,
        94228,  400468, 54048,  66816,  131304, 61176,  109776, 54528,
        129120, 205780, 94804,  92500,  79584,  181928, 208372};

/**
 * @brief One solid's edge geometry, held in the composed runtime's
 *        persistent state for the compiled Tiling source.
 * @details The block sits at a fixed position inside the runtime state, so a
 * rebuild reuses it whatever the runtime allocated after it. It holds the
 * narrow index alone while the current and previous frames both draw thin
 * bands, so every column-time band between them reads it, and the coarse
 * index otherwise.
 */
class SphereTilingEdgeStore {
public:
  /// Geometry block capacity.
  static constexpr size_t BYTES = 164 * 1024;
  /// Block bytes left for device rounding to admit a few more index candidates.
  static constexpr size_t ROUNDING_SLACK = 2 * 1024;

  /** @brief Whether a solid's edge geometry fits the block.
   * @param solid Registry solid ordinal.
   * @return True for a registry solid within BYTES less ROUNDING_SLACK. */
  static constexpr bool admits(unsigned solid) {
    return solid < SPHERE_TILING_EDGE_BYTES.size() &&
           SPHERE_TILING_EDGE_BYTES[solid] + ROUNDING_SLACK <= BYTES;
  }

  /** @brief Whether two parameter sets select the same geometry.
   * @param a First tiling parameters.
   * @param b Second tiling parameters.
   * @return True when both name the same solid. */
  static constexpr bool same_geometry(const SphericalTilingSourceParams &a,
                                      const SphericalTilingSourceParams &b) {
    return a.solid == b.solid;
  }

  /** @brief Builds the parameters' solid unless it is already built.
   * @param params Live tiling parameters; the solid must be admitted.
   * @return True when this call rebuilt the geometry. */
  HS_FLASH_MEMBER bool refresh(const SphericalTilingSourceParams &params) {
    const bool thin = thin_band(params);
    const bool narrow = thin && (previous_thin || built == NO_SOLID);
    previous_thin = thin;
    if (params.solid == built && narrow == narrow_built)
      return false;
    // no death case: valid() refuses unadmitted solids before they go live.
    HS_CHECK(admits(params.solid),
             "sphere tiling: solid %u exceeds the edge store",
             static_cast<unsigned>(params.solid));
    built = NO_SOLID;
    using Indices = SphereTilingGeometry::Indices;
    geometry = SphereTilingGeometry::build(
        params.solid, block, scratch_arena_a, scratch_arena_b,
        narrow ? Indices::NARROW : Indices::COARSE);
    built = params.solid;
    narrow_built = narrow;
    return true;
  }

  /** @brief Whether parameters draw a band the narrow index serves.
   * @details Its margin keeps every band interpolated between two such
   * parameter sets inside the narrow limit: the outer angle is linear in the
   * thickness and convex in the log-interpolated softness.
   * @param params Tiling parameters.
   * @return True when the band's outer angle is within the narrow limit. */
  static bool thin_band(const SphericalTilingSourceParams &params) {
    return SphereTilingEdgeBand(params.line_thickness, params.line_softness)
               .outer <=
           SphereTilingNarrowIndex::OUTER * (1.0f - THIN_BAND_MARGIN);
  }

  /** @brief The built solid's edges, which every frame samples.
   * @return Indexed edges, or an empty view before the first build. */
  SphereTilingEdgeView edges() const {
    return built == NO_SOLID ? SphereTilingEdgeView{} : geometry.edge_view();
  }

  /** @brief Registry ordinal of the built solid.
   * @return The solid, or UINT8_MAX before the first build. */
  uint8_t built_solid() const { return built; }

  /** @brief Points a frame at this store.
   * @tparam Frame Composed frame state.
   * @param frame Frame to bind. */
  template <typename Frame> void bind(Frame &frame) const {
    frame.sphere_tiling = this;
  }

  /** @brief Bytes the built geometry occupies.
   * @return Used block bytes, or zero before the first build. */
  size_t bytes() const { return built == NO_SOLID ? 0 : geometry.bytes(); }

private:
  static constexpr uint8_t NO_SOLID = UINT8_MAX;
  static constexpr float THIN_BAND_MARGIN = 1.0e-4f;

  alignas(std::max_align_t) std::array<uint8_t, BYTES> block;
  SphereTilingGeometry geometry;
  uint8_t built = NO_SOLID;
  bool narrow_built = false;
  bool previous_thin = false;
};

static_assert(Solids::NUM_ENTRIES < UINT8_MAX);

/** @brief Checks parameters for the compiled edge-only source.
 * @param params Candidate tiling parameters.
 * @return True if the fields are valid, edge output is selected, and the
 * solid's geometry fits the edge store. */
inline bool valid(const SphericalTilingSourceParams &params) {
  return Fields::valid(params) &&
         params.output == static_cast<uint8_t>(TilingOutput::EDGES) &&
         SphereTilingEdgeStore::admits(params.solid);
}

/** @brief Samples the nearest edge of the runtime's stored spherical tiling.
 * @tparam State Provider of tiling parameters.
 */
template <typename State> struct StoredTilingEdges : ApproximationDefaults {
  using Binding =
      typename State::Binding; ///< Pipeline binding supplied by the provider.
  using FrameState =
      typename State::FrameState; ///< Frame state read by the provider.

  static_assert(
      !requires(const FrameState &frame) { frame.source_morph; },
      "the stored tiling holds one endpoint and cannot morph");

  /// Edge geometry and coverage controls for one frame.
  struct Prepared {
    SphereTilingEdgeView edges; ///< View into the runtime's edge store.
    SphereTilingEdgeBand band{0.0f,
                              0.05f}; ///< Angular stroke and fringe limits.
  };

  /** @brief Resolves the stored geometry and the frame's band.
   * @details Columns sampled before a solid change in their frame draw the
   * newly built solid.
   * @param frame Current parameter state and edge store.
   * @return The stored edges, or none without a store. */
  static Prepared prepare(const FrameState &frame) {
    const auto &params = State::params(frame);
    return {frame.sphere_tiling ? frame.sphere_tiling->edges()
                                : SphereTilingEdgeView{},
            SphereTilingEdgeBand(params.line_thickness, params.line_softness)};
  }

  /** @brief Samples the nearest stored edges in the input direction.
   * @param input Sphere direction and carried path length.
   * @param source Prepared edge geometry and band.
   * @return Arc fraction and coverage, with the input direction and path
   * length. */
  static HS_HOT_INLINE FieldSample sample(const SphereSample &input,
                                          const FrameState &,
                                          const Prepared &source) {
    const auto direction = Pattern::SphereDirection::admit(input.dir);
    const auto projection = math::project_cube_face(direction.value());
    const auto edge =
        sample_edges(direction.value(), source, projection,
                     SphereTilingNarrowIndex::cell_for(projection));
    return {edge.value, edge.coverage, input.dir, input.path_length};
  }

private:
  static HS_HOT_ITCM_MEMBER SphereTilingEdgeValue
  sample_edges(const math::Vector &direction, const Prepared &source,
               const math::CubeFaceProjection &projection, int narrow_cell) {
    return source.edges.sample(direction, source.band, projection, narrow_cell);
  }
};

} // namespace Source

/** @brief Stored tiling edge source. @tparam B Frame binding. */
template <typename B>
struct SourcePolicyFor<Source::SphericalTilingSourceParams, B> {
  using Type = Source::StoredTilingEdges<
      SourceProvider<B, Source::SphericalTilingSourceParams>>;
};

/** @brief The runtime keeps the Tiling source's edge geometry. */
template <> struct SourceStoreFor<Source::SphericalTilingSourceParams> {
  using Type = Source::SphereTilingEdgeStore; ///< Persistent edge store.
};

} // namespace Pullback
