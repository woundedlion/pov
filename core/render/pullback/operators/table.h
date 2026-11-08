/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "core/render/pullback/operators.h"
#include "core/render/pullback/operators/snapshot.h"

/**
 * @file table.h
 * @brief The chain interpreter's operator table; the catalog is pinned to it.
 */

namespace Pullback {

namespace Interp {

/** @brief Every chain operator model, in catalog order. */
using OperatorModels = ::Pullback::Detail::TypeList<
    Op::Rotate, Op::DisplaceCurl, Op::DisplaceDirect, Op::DisplaceRipple,
    Op::DisplacePotential, Op::DisplaceRossbyHaurwitz, Op::DisplaceSymmetric,
    Op::LensGlitch, Op::LensTwist, Op::LensMobius, Op::LensKaleidoscope,
    Op::ProjectStereographic, Op::ProjectFoldedSinusoidal,
    Op::ProjectEquirectangular, Op::ProjectGnomonic, Op::ProjectPeirceV3,
    Op::ProjectPeirceSquareFastV3, Op::ProjectBonneV3, Op::ProjectAiroceanV3,
    Op::ProjectTetrahedral, Op::ProjectCubic, Op::ProjectOctahedral,
    Op::ProjectDodecahedral, Op::ProjectIcosahedral, Op::WarpAffineV4,
    Op::WarpWaveShear, Op::WarpVortex, Op::WarpVectorNoise, Op::WarpMirrorTile,
    Op::WarpPolarChart, Op::WarpCurlFlow, Op::SampleGridV4,
    Op::SampleTwinWaveV3, Op::SampleBands, Op::SampleSpiralV3,
    Op::SampleTessellation, Op::SampleProjectedNoise, Op::SampleSphericalNoise,
    Op::SampleSphericalCells, Op::TransferRidge, Op::TransferIsoContour,
    Op::TransferSmoothBands, Op::CoverageValueCutout,
    Op::ColorizeGeneratedPaletteV4, Op::ModulateColorV1>;

/** @brief One descriptor per model of @p Models.
 *  @tparam Models Operator models.
 *  @return The models' descriptors, in list order. */
template <typename... Models>
consteval std::array<OperatorDescriptor, sizeof...(Models)>
make_operator_table(::Pullback::Detail::TypeList<Models...>) {
  return {make_operator_descriptor<Models>()...};
}

/** @brief Every chain operator's descriptor, one per operator model. */
inline constexpr auto OPERATOR_TABLE = make_operator_table(OperatorModels{});

/**
 * @brief Whether every table operator id is unique.
 * @return False on a duplicate id.
 */
consteval bool operator_ids_unique() {
  for (size_t i = 0; i < OPERATOR_TABLE.size(); ++i)
    for (size_t j = 0; j < i; ++j)
      if (std::string_view(OPERATOR_TABLE[i].operator_id) ==
          OPERATOR_TABLE[j].operator_id)
        return false;
  return true;
}

/**
 * @brief Whether every table display name is unique.
 * @return False on a duplicate name.
 */
consteval bool operator_names_unique() {
  for (size_t i = 0; i < OPERATOR_TABLE.size(); ++i)
    for (size_t j = 0; j < i; ++j)
      if (std::string_view(OPERATOR_TABLE[i].display_name) ==
          OPERATOR_TABLE[j].display_name)
        return false;
  return true;
}

static_assert(operator_ids_unique(),
              "chain operator table: duplicate operator id");
static_assert(operator_names_unique(),
              "chain operator table: duplicate operator name");

/** @brief The table entry named @p operator_id, or null.
 *  @param operator_id Operator id to look up.
 *  @return The descriptor, or null when absent. */
inline const OperatorDescriptor *find_operator(std::string_view operator_id) {
  for (const OperatorDescriptor &op : OPERATOR_TABLE)
    if (operator_id == op.operator_id)
      return &op;
  return nullptr;
}

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
