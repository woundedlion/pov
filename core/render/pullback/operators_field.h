/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "render/pullback/material.h"
#include "render/pullback/operators_common.h"

/**
 * @file operators_field.h
 * @brief FIELD-endomorphism operator models: the value transfers and the
 *        value-dependent coverage cut.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Shared shape of the stateless FIELD endomorphisms. */
template <typename Derived, typename ParamsT>
struct FieldEndoModel : StatelessModel {
  using Input = FieldSample;
  using Output = FieldSample;
  using Params = ParamsT;
  struct Prepared {};

  static Prepared prepare(const FrameContext &, const Params &,
                          const StatelessModel::State &) {
    return {};
  }
};

/** @brief FIELD endomorphism: the unit-bell ridge transfer. */
struct TransferRidge : FieldEndoModel<TransferRidge, Transfer::NoValueParams> {
  static constexpr const char *ID = "field.transfer.ridge.v2";
  static constexpr const char *NAME = "Ridge Transfer";

  static FieldSample run(const FieldSample &input, const FrameContext &,
                         const Params &, const Prepared &) {
    return Kernel::transfer(input, unit_bell(input.value));
  }
};

/** @brief Parameter family of field.transfer.iso-contour.v2. */
struct IsoContourChainParams : Transfer::IsoValueParams {
  static constexpr auto FIELDS = concat_fields<IsoContourChainParams>(
      Transfer::IsoValueParams::FIELDS,
      std::array<Field<IsoContourChainParams>, 0>{});
};
static_assert(field_ids_unique<IsoContourChainParams>());

/** @brief FIELD endomorphism: the iso-band transfer. */
struct TransferIsoContour
    : FieldEndoModel<TransferIsoContour, IsoContourChainParams> {
  static constexpr const char *ID = "field.transfer.iso-contour.v2";
  static constexpr const char *NAME = "Iso Contour";

  static FieldSample run(const FieldSample &input, const FrameContext &,
                         const Params &params, const Prepared &) {
    return Kernel::transfer(
        input,
        Transfer::iso_contour(input.value, params.iso_level, params.iso_width));
  }
};

/** @brief Parameter family of field.transfer.smooth-bands.v2. */
struct SmoothBandsChainParams {
  float band_count = 4.0f; /**< Cosine bands over the unit value. */
  float band_phase = 0.0f; /**< Band offset, in radians. */

  static constexpr auto FIELDS = std::array{
      Field<SmoothBandsChainParams>{
          "band-count", &SmoothBandsChainParams::band_count, "Band Count", 1.0f,
          32.0f, FieldCurve::SNAP},
      Field<SmoothBandsChainParams>{
          "band-phase", &SmoothBandsChainParams::band_phase, "Band Phase", 0.0f,
          TWO_PI_F, FieldCurve::SHORTEST_PERIODIC},
  };
};
static_assert(field_ids_unique<SmoothBandsChainParams>());

/** @brief FIELD endomorphism: the cosine banding transfer. */
struct TransferSmoothBands
    : FieldEndoModel<TransferSmoothBands, SmoothBandsChainParams> {
  static constexpr const char *ID = "field.transfer.smooth-bands.v2";
  static constexpr const char *NAME = "Smooth Bands";

  static FieldSample run(const FieldSample &input, const FrameContext &,
                         const Params &params, const Prepared &) {
    return Kernel::transfer(input, Transfer::smooth_bands(input.value,
                                                          params.band_count,
                                                          params.band_phase));
  }
};

/** @brief Parameter family of field.coverage.value-cutout.v2. */
struct ValueCutoutChainParams : ValueCoverage::CutoutValueParams {
  static constexpr auto FIELDS = concat_fields<ValueCutoutChainParams>(
      ValueCoverage::CutoutValueParams::FIELDS,
      std::array<Field<ValueCutoutChainParams>, 0>{});
};
static_assert(field_ids_unique<ValueCutoutChainParams>());

/** @brief FIELD endomorphism: the value-dependent coverage cut. */
struct CoverageValueCutout
    : FieldEndoModel<CoverageValueCutout, ValueCutoutChainParams> {
  static constexpr const char *ID = "field.coverage.value-cutout.v2";
  static constexpr const char *NAME = "Value Cutout";

  static FieldSample run(const FieldSample &input, const FrameContext &,
                         const Params &params, const Prepared &) {
    return Kernel::coverage(
        input, ValueCoverage::value_cutout(input.value, params.cutout_threshold,
                                           params.cutout_softness));
  }
};

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
