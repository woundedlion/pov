/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "render/pullback/operators.h"

/**
 * @file operator_table.h
 * @brief The chain interpreter's operator table: the C++ ground truth the
 *        catalog is pinned to and setShaderChain resolves against.
 */

namespace Pullback {

namespace Interp {

inline constexpr std::array<OperatorDescriptor, 34> OPERATOR_TABLE{
    make_operator_descriptor<Op::Rotate>(),
    make_operator_descriptor<Op::DisplaceCurl>(),
    make_operator_descriptor<Op::DisplaceDirect>(),
    make_operator_descriptor<Op::LensGlitch>(),
    make_operator_descriptor<Op::LensTwist>(),
    make_operator_descriptor<Op::LensMobius>(),
    make_operator_descriptor<Op::LensKaleidoscope>(),
    make_operator_descriptor<Op::ProjectStereographic>(),
    make_operator_descriptor<Op::ProjectFoldedSinusoidal>(),
    make_operator_descriptor<Op::ProjectEquirectangular>(),
    make_operator_descriptor<Op::ProjectGnomonic>(),
    make_operator_descriptor<Op::ProjectPeirce>(),
    make_operator_descriptor<Op::ProjectPeirceSquareFast>(),
    make_operator_descriptor<Op::ProjectBonne>(),
    make_operator_descriptor<Op::ProjectAirocean>(),
    make_operator_descriptor<Op::WarpAffine>(),
    make_operator_descriptor<Op::WarpWaveShear>(),
    make_operator_descriptor<Op::WarpVectorNoise>(),
    make_operator_descriptor<Op::WarpMirrorTile>(),
    make_operator_descriptor<Op::WarpPolarChart>(),
    make_operator_descriptor<Op::WarpCurlFlow>(),
    make_operator_descriptor<Op::SampleGrid>(),
    make_operator_descriptor<Op::SampleTwinWave>(),
    make_operator_descriptor<Op::SampleRings>(),
    make_operator_descriptor<Op::SampleSpiral>(),
    make_operator_descriptor<Op::SampleLattice>(),
    make_operator_descriptor<Op::SampleProjectedNoise>(),
    make_operator_descriptor<Op::SampleSphericalNoise>(),
    make_operator_descriptor<Op::TransferLinear>(),
    make_operator_descriptor<Op::TransferRidge>(),
    make_operator_descriptor<Op::TransferIsoContour>(),
    make_operator_descriptor<Op::TransferSmoothBands>(),
    make_operator_descriptor<Op::CoverageValueCutout>(),
    make_operator_descriptor<Op::ColorizeGeneratedPalette>(),
};

consteval bool operator_ids_unique() {
  for (size_t i = 0; i < OPERATOR_TABLE.size(); ++i)
    for (size_t j = 0; j < i; ++j)
      if (std::string_view(OPERATOR_TABLE[i].operator_id) ==
          OPERATOR_TABLE[j].operator_id)
        return false;
  return true;
}

/** Per-op monotonicity is a table invariant: adjacency over monotone
    operators yields a monotone chain, so compile() never re-walks it. */
consteval bool operator_table_monotone() {
  for (const OperatorDescriptor &op : OPERATOR_TABLE)
    if (static_cast<uint8_t>(op.input) > static_cast<uint8_t>(op.output))
      return false;
  return true;
}

static_assert(operator_ids_unique(),
              "chain operator table: duplicate operator id");
static_assert(operator_table_monotone(),
              "chain operator table: an operator may not decrease its family "
              "rank");

/** @brief The table entry named @p operator_id, or null. */
inline const OperatorDescriptor *find_operator(std::string_view operator_id) {
  for (const OperatorDescriptor &op : OPERATOR_TABLE)
    if (operator_id == op.operator_id)
      return &op;
  return nullptr;
}

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
