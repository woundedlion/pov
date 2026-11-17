/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/math/scalar_kernels.h"
#include "core/math/3dmath.h"
#include "core/math/scalar_noise.h"
#include "core/platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "core/render/pullback/material.h"
#include "core/render/pullback/operators/common.h"

/**
 * @file field.h
 * @brief FIELD transfer and deterministic or stochastic coverage operators.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Shared shape of the stateless FIELD endomorphisms. */
template <typename ParamsT> struct FieldEndoModel : StatelessModel {
  using Input = FieldSample;  ///< FIELD in.
  using Output = FieldSample; ///< FIELD out.
  using Params = ParamsT;     ///< The transfer's param family.
  struct Prepared {};

  /** @brief Empty prepared block. @return `{}`. */
  static Prepared prepare(const FrameContext &, const Params &,
                          const StatelessModel::State &) {
    return {};
  }
};

/** @brief FIELD endomorphism: the unit-bell ridge transfer. */
struct TransferRidge : FieldEndoModel<Transfer::NoValueParams> {
  static constexpr const char *ID = "field.transfer.ridge.v2";
  static constexpr const char *NAME = "Ridge Transfer";

  static FieldSample run(const FieldSample &input, const FrameContext &,
                         const Params &, const Prepared &) {
    return Kernel::transfer(input, math::unit_bell(input.value));
  }
};

/** @brief Parameter family of field.transfer.iso-contour.v2. */
using IsoContourChainParams = Transfer::IsoValueParams;

/** @brief FIELD endomorphism: the iso-band transfer. */
struct TransferIsoContour : FieldEndoModel<IsoContourChainParams> {
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
using SmoothBandsChainParams = Transfer::SmoothBandsParams;

/** @brief FIELD endomorphism: the cosine banding transfer. */
struct TransferSmoothBands : FieldEndoModel<SmoothBandsChainParams> {
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
using ValueCutoutChainParams = ValueCoverage::CutoutValueParams;

/** @brief FIELD endomorphism: the value-dependent coverage cut. */
struct CoverageValueCutout : FieldEndoModel<ValueCutoutChainParams> {
  static constexpr const char *ID = "field.coverage.value-cutout.v2";
  static constexpr const char *NAME = "Value Cutout";

  static FieldSample run(const FieldSample &input, const FrameContext &,
                         const Params &params, const Prepared &) {
    return Kernel::coverage(
        input, ValueCoverage::value_cutout(input.value, params.cutout_threshold,
                                           params.cutout_softness));
  }
};

/** @brief Probability response distribution. */
enum class ProbabilityDistribution : uint8_t {
  UNIFORM, ///< Uniform threshold.
  POISSON, ///< Poisson events.
  NORMAL   ///< Normal threshold.
};

/** @brief Parameters of field.coverage.probability.v1. */
struct ProbabilityCoverageParams {
  float field_min =
      0.0f; ///< Drive endpoint; the lower field endpoint maps to zero.
  float field_max =
      1.0f; ///< Drive endpoint; the upper field endpoint maps to one.
  float probability_min =
      0.0f; ///< Probability endpoint; the pair is sorted into bounds.
  float probability_max =
      1.0f;          ///< Probability endpoint; the pair is sorted into bounds.
  float rate = 1.0f; ///< Poisson events per frame at unit drive.
  float mean = 0.5f; ///< Normal threshold mean in drive units.
  float deviation = 0.2f;   ///< Normal threshold standard deviation.
  float seed = 0.0f;        ///< Additional integer random stream seed.
  uint8_t distribution = 0; ///< Uniform, Poisson or normal response.

  /** @brief Distribution selector spellings. */
  static constexpr const char *DISTRIBUTIONS[] = {"uniform", "poisson",
                                                  "normal"};
  static constexpr auto FIELDS = std::to_array<
      Field<ProbabilityCoverageParams>>({
      {"field-min", &ProbabilityCoverageParams::field_min, "Field Min", 0, 1},
      {"field-max", &ProbabilityCoverageParams::field_max, "Field Max", 0, 1},
      {"probability-min", &ProbabilityCoverageParams::probability_min,
       "Probability Min", 0, 1},
      {"probability-max", &ProbabilityCoverageParams::probability_max,
       "Probability Max", 0, 1},
      {"rate",
       &ProbabilityCoverageParams::rate,
       "Poisson Rate",
       0,
       16,
       FieldCurve::LERP,
       FieldGate::ALWAYS,
       {"distribution",
        live_values(static_cast<uint8_t>(ProbabilityDistribution::POISSON))}},
      {"mean",
       &ProbabilityCoverageParams::mean,
       "Normal Mean",
       0,
       1,
       FieldCurve::LERP,
       FieldGate::ALWAYS,
       {"distribution",
        live_values(static_cast<uint8_t>(ProbabilityDistribution::NORMAL))}},
      {"deviation",
       &ProbabilityCoverageParams::deviation,
       "Normal Deviation",
       0.01f,
       1,
       FieldCurve::LOG_POSITIVE,
       FieldGate::ALWAYS,
       {"distribution",
        live_values(static_cast<uint8_t>(ProbabilityDistribution::NORMAL))}},
      {"seed", &ProbabilityCoverageParams::seed, "Seed", 0, 65535,
       FieldCurve::SNAP},
  });
  static constexpr auto TOPOLOGY =
      std::to_array<TopologyField<ProbabilityCoverageParams>>({
          {"distribution", &ProbabilityCoverageParams::distribution,
           DISTRIBUTIONS, 0},
      });
};
static_assert(std::size(ProbabilityCoverageParams::DISTRIBUTIONS) ==
              static_cast<size_t>(ProbabilityDistribution::NORMAL) + 1);
static_assert(field_ids_unique<ProbabilityCoverageParams>());
static_assert(field_defaults_in_range<ProbabilityCoverageParams>());

/** @brief Instance seed and discrete frame clock of probability coverage. */
struct ProbabilityClockState {
  uint32_t seed = 0;  ///< Instance random stream.
  uint32_t frame = 0; ///< Frames advanced, modulo 2^32.
};

/** @brief Empty interval for a discrete frame clock.
 * @return An empty interval. */
inline EmptyState temporal_interval(const ProbabilityClockState &,
                                    const ProbabilityClockState &) {
  return {};
}

/** @brief Holds the current frame's random stream across column times. */
inline void temporal_sample(ProbabilityClockState &, const EmptyState &,
                            const ProbabilityCoverageParams &, float) {}

/** @brief FIELD coverage: independent Bernoulli trials driven by field value. */
struct CoverageProbability : ValueStateModel<ProbabilityClockState> {
  static constexpr const char *ID = "field.coverage.probability.v1";
  static constexpr const char *NAME = "Probability Coverage";
  using Input = FieldSample;
  using Output = FieldSample;
  using Params = ProbabilityCoverageParams;
  /** @brief Frame seed and normal-distribution constants. */
  struct Prepared {
    uint32_t seed;      ///< Instance and authored stream seed.
    uint32_t frame;     ///< Discrete frame index.
    float normal_low;   ///< Error function at zero drive.
    float normal_span;  ///< Error-function range across unit drive.
    float normal_scale; ///< 1 / (deviation * sqrt(2)).
  };

  static void init(State &state, InstanceId id) { state.seed = id.stable_hash; }
  static void advance(State &state, const Params &) { ++state.frame; }

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    Prepared out{state.seed ^ static_cast<uint32_t>(params.seed), state.frame,
                 0, 1, 1};
    if (params.distribution ==
        static_cast<uint8_t>(ProbabilityDistribution::NORMAL)) {
      out.normal_scale = 1.0f / (params.deviation * sqrtf(2.0f));
      out.normal_low = erff(-params.mean * out.normal_scale);
      out.normal_span =
          erff((1.0f - params.mean) * out.normal_scale) - out.normal_low;
    }
    return out;
  }

  /** @brief Lighting probability for a normalized field sample.
   * @param value Input field value.
   * @param params Distribution and range controls.
   * @param prepared Frame constants.
   * @return Probability in [0, 1], nondecreasing in @p value. */
  static float probability(float value, const Params &params,
                           const Prepared &prepared) {
    const float low = std::min(params.field_min, params.field_max);
    const float high = std::max(params.field_min, params.field_max);
    float drive = high > low
                      ? std::clamp((value - low) / (high - low), 0.0f, 1.0f)
                      : float(value >= high);
    if (params.distribution ==
        static_cast<uint8_t>(ProbabilityDistribution::POISSON))
      drive = -expm1f(-params.rate * drive);
    else if (params.distribution ==
             static_cast<uint8_t>(ProbabilityDistribution::NORMAL))
      drive = std::clamp((erff((drive - params.mean) * prepared.normal_scale) -
                          prepared.normal_low) /
                             prepared.normal_span,
                         0.0f, 1.0f);
    const float minimum =
        std::min(params.probability_min, params.probability_max);
    const float maximum =
        std::max(params.probability_min, params.probability_max);
    return minimum + (maximum - minimum) * drive;
  }

  static FieldSample run(const FieldSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    uint32_t point = hs::FNV1A_32_BASIS;
    for (const float coordinate :
         {input.sphere.x, input.sphere.y, input.sphere.z}) {
      const auto cell = static_cast<uint32_t>(
          (std::clamp(coordinate, -1.0f, 1.0f) + 1.0f) * 32768.0f);
      point = (point ^ cell) * hs::FNV1A_32_PRIME;
    }
    const float draw =
        math::hash01(point, prepared.seed + prepared.frame * 0x9e3779b9U);
    return Kernel::coverage(
        input, float(draw < probability(input.value, params, prepared)));
  }
};

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
