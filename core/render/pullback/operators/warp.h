/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "render/pullback/operators/common.h"
#include "render/pullback/warp.h"

/**
 * @file warp.h
 * @brief PLANE-endomorphism operator models over the shared warp kernels.
 */

namespace Pullback {

namespace Interp {

namespace Op {

using WarpEnvelope = Warp::Envelope;

inline constexpr auto &WARP_ENVELOPE_IDS = Warp::ENVELOPE_IDS;

/**
 * @brief Bounds a warp operator's envelope enum8.
 * @details warp_envelope() relies on this check and carries no per-pixel guard.
 */
inline void check_warp_envelope(uint8_t envelope) {
  HS_CHECK(envelope <= static_cast<uint8_t>(WarpEnvelope::EDGE_FADE),
           "warp operator: invalid envelope");
}

/** @brief The envelope switch over the shared warp envelope kernel. */
inline float warp_envelope(uint8_t envelope,
                           const ProjectionProvenance &provenance,
                           float edge_width) {
  return Warp::envelope(provenance, edge_width,
                        static_cast<WarpEnvelope>(envelope));
}

/** @brief Phase clock of the single-clock warp operators. */
struct WarpPhaseState {
  float phase = 0.0f;
};

/** @brief Parameter family of warp.affine.v3; translations are lattice cells. */
struct AffineWarpParams : Warp::AffineParams {
  float lattice_period = 1.0f;
  static constexpr auto FIELDS = concat_fields<AffineWarpParams>(
      Warp::AffineParams::FIELDS,
      std::array{Field<AffineWarpParams>{
          "lattice-period", &AffineWarpParams::lattice_period, "Lattice Period",
          1.0f / 64.0f, 100.0f, FieldCurve::LOG_POSITIVE}});
};
static_assert(field_ids_unique<AffineWarpParams>());
static_assert(
    appended_block_size_matches<AffineWarpParams, Warp::AffineParams, 4>(),
    "appended parameter block must have the expected rounded size");
static_assert(field_defaults_in_range<AffineWarpParams>());

/** @brief Phase clock plus the accumulated frame rotation of warp.affine.v3. */
struct AffineClockState {
  float phase = 0.0f;
  float rotation = 0.0f;
};

/** @brief PLANE endomorphism: the oscillating affine frame change. */
struct WarpAffineV3 : ValueStateModel<AffineClockState> {
  static constexpr const char *ID = "warp.affine.v3";
  static constexpr const char *NAME = "Affine Warp Extended Period";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = AffineWarpParams;
  using Prepared = Warp::PreparedAffineSlot;

  static void advance(State &state, const Params &params) {
    state.phase = math::wrap_t(state.phase + params.speed);
    Warp::advance_affine_rotation(state.rotation, params);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    return Warp::prepare(params, state.phase, state.rotation,
                         params.lattice_period);
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &, const Prepared &prepared) {
    return Kernel::warp(input,
                        Warp::affine_frame(input.coords, prepared, true));
  }
};

/** @brief Parameter family of warp.wave-shear.v2. */
struct WaveShearWarpParams : Warp::WaveShearParams {
  uint8_t envelope = static_cast<uint8_t>(WarpEnvelope::FLAT);

  static constexpr auto TOPOLOGY = std::array{
      TopologyField<WaveShearWarpParams>{
          "envelope", &WaveShearWarpParams::envelope, WARP_ENVELOPE_IDS,
          static_cast<uint8_t>(WarpEnvelope::FLAT)},
  };
};
static_assert(field_ids_unique<WaveShearWarpParams>());
static_assert(appended_block_size_matches<WaveShearWarpParams,
                                          Warp::WaveShearParams, 1>(),
              "appended parameter block must have the expected rounded size");
static_assert(field_defaults_in_range<WaveShearWarpParams>());

/** @brief The wave shear's prepared block: the field frame plus the phase the
    kernel consumes per sample. */
struct PreparedWaveShear {
  Warp::PreparedRotation rotation;
  float phase;
};

/** @brief PLANE endomorphism: the travelling sine shear. */
struct WarpWaveShear : PhaseClockModel<WarpPhaseState> {
  static constexpr const char *ID = "warp.wave-shear.v2";
  static constexpr const char *NAME = "Wave Shear";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = WaveShearWarpParams;
  using Prepared = PreparedWaveShear;

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_warp_envelope(params.envelope);
    return {Warp::prepare(params, state.phase), state.phase};
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    const float amplitude =
        params.strength *
        warp_envelope(params.envelope, input.provenance, params.edge_width);
    return Kernel::warp(input,
                        Warp::wave_shear(input.coords, params, prepared.phase,
                                         amplitude, prepared.rotation, true));
  }
};

/** @brief Parameter family of warp.vortex.v2. */
using VortexWarpParams = Warp::VortexParams;

/** @brief PLANE endomorphism: the orbiting radial vortex. */
struct WarpVortex : PhaseClockModel<WarpPhaseState> {
  static constexpr const char *ID = "warp.vortex.v2";
  static constexpr const char *NAME = "Vortex";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = VortexWarpParams;
  using Prepared = Warp::PreparedVortexSlot;

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    return Warp::prepare(params, state.phase);
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &, const Prepared &prepared) {
    return Kernel::warp(input, Warp::vortex(input.coords, prepared, true));
  }
};

/** @brief Parameter family of warp.vector-noise.v2. */
struct VectorNoiseWarpParams : Warp::VectorNoiseParams {
  uint8_t basis = static_cast<uint8_t>(math::NoiseBasis::SIMPLEX);
  uint8_t envelope = static_cast<uint8_t>(WarpEnvelope::FLAT);

  static constexpr auto TOPOLOGY = std::array{
      TopologyField<VectorNoiseWarpParams>{
          "basis", &VectorNoiseWarpParams::basis, NOISE_BASIS_IDS,
          static_cast<uint8_t>(math::NoiseBasis::SIMPLEX)},
      TopologyField<VectorNoiseWarpParams>{
          "envelope", &VectorNoiseWarpParams::envelope, WARP_ENVELOPE_IDS,
          static_cast<uint8_t>(WarpEnvelope::FLAT)},
  };
};
static_assert(field_ids_unique<VectorNoiseWarpParams>());
static_assert(appended_block_size_matches<VectorNoiseWarpParams,
                                          Warp::VectorNoiseParams, 2>(),
              "appended parameter block must have the expected rounded size");
static_assert(field_defaults_in_range<VectorNoiseWarpParams>());

/** @brief The vector-noise warp's prepared block: the owned noise field plus
    the vector frame and loop point. */
struct PreparedVectorNoiseWarp {
  const FastNoiseLite *noise;
  Warp::PreparedVectorNoiseSlot slot;
};

/** @brief PLANE endomorphism: the noise-vector displacement. */
struct WarpVectorNoise : PhaseClockModel<NoisePhaseState> {
  static constexpr const char *ID = "warp.vector-noise.v2";
  static constexpr const char *NAME = "Vector Noise";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = VectorNoiseWarpParams;
  using Prepared = PreparedVectorNoiseWarp;

  static void init(State &state, InstanceId id) { init_noise_phase(state, id); }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_warp_envelope(params.envelope);
    check_noise_basis(params.basis);
    return {&state.noise, Warp::prepare(params, state.phase)};
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    if (!noise_plane_in_domain(input.coords))
      return input;
    const float amplitude =
        params.strength *
        warp_envelope(params.envelope, input.provenance, params.edge_width);
    return Kernel::warp(
        input,
        Warp::vector_noise(input.coords, params, amplitude, *prepared.noise,
                           static_cast<math::NoiseBasis>(params.basis),
                           prepared.slot, true));
  }
};

/** @brief Parameter family of warp.mirror-tile.v2. */
using MirrorWarpParams = Warp::MirrorParams;

/** @brief PLANE endomorphism: the mirrored tiling fold. */
struct WarpMirrorTile : PhaseClockModel<WarpPhaseState> {
  static constexpr const char *ID = "warp.mirror-tile.v2";
  static constexpr const char *NAME = "Mirror Tile";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = MirrorWarpParams;
  using Prepared = Warp::PreparedMirrorSlot;

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    return Warp::prepare(params, state.phase);
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    return Kernel::warp(
        input, Warp::mirror_tile(input.coords, params, prepared, true));
  }
};

enum class PolarMode : uint8_t { LINEAR = 0, LOGARITHMIC = 1 };

inline constexpr const char *POLAR_MODE_IDS[] = {"linear", "logarithmic"};
static_assert(std::size(POLAR_MODE_IDS) ==
              static_cast<size_t>(PolarMode::LOGARITHMIC) + 1);

/** @brief Harmonic topology values h1..h16, indexed by harmonic - 1. */
inline constexpr const char *POLAR_HARMONIC_IDS[] = {
    "h1", "h2",  "h3",  "h4",  "h5",  "h6",  "h7",  "h8",
    "h9", "h10", "h11", "h12", "h13", "h14", "h15", "h16"};
static_assert(std::size(POLAR_HARMONIC_IDS) == Warp::MAX_POLAR_HARMONIC);

/** @brief Parameter family of warp.polar-chart.v2. */
struct PolarChartParams : Warp::PolarParams {
  uint8_t mode = static_cast<uint8_t>(PolarMode::LINEAR);
  uint8_t harmonic = 0; /**< Harmonic value index; harmonic = index + 1. */

  static constexpr auto TOPOLOGY = std::array{
      TopologyField<PolarChartParams>{"mode", &PolarChartParams::mode,
                                      POLAR_MODE_IDS,
                                      static_cast<uint8_t>(PolarMode::LINEAR)},
      TopologyField<PolarChartParams>{"harmonic", &PolarChartParams::harmonic,
                                      POLAR_HARMONIC_IDS, 0},
  };
};
static_assert(field_ids_unique<PolarChartParams>());
static_assert(
    appended_block_size_matches<PolarChartParams, Warp::PolarParams, 2>(),
    "appended parameter block must have the expected rounded size");
static_assert(field_defaults_in_range<PolarChartParams>());

/** @brief PLANE endomorphism: the polar chart change. */
struct WarpPolarChart : PhaseClockModel<WarpPhaseState> {
  static constexpr const char *ID = "warp.polar-chart.v2";
  static constexpr const char *NAME = "Polar Chart";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = PolarChartParams;
  struct Prepared {
    float phase;
  };

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    HS_CHECK(params.mode <= static_cast<uint8_t>(PolarMode::LOGARITHMIC),
             "warp.polar-chart: invalid polar mode");
    HS_CHECK(params.harmonic < Warp::MAX_POLAR_HARMONIC,
             "warp.polar-chart: invalid harmonic");
    return {state.phase};
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    return Kernel::warp(
        input, Warp::polar_chart(input.coords, params, prepared.phase,
                                 static_cast<PolarMode>(params.mode) ==
                                     PolarMode::LOGARITHMIC,
                                 static_cast<uint8_t>(params.harmonic + 1)));
  }
};

enum class CurlIntegrator : uint8_t { EULER1, MIDPOINT2, MIDPOINT4 };

inline constexpr uint8_t curl_intervals(CurlIntegrator integrator) {
  static_assert(Warp::Euler1::INTERVALS ==
                (1U << static_cast<uint8_t>(CurlIntegrator::EULER1)));
  static_assert(Warp::Midpoint2::INTERVALS ==
                (1U << static_cast<uint8_t>(CurlIntegrator::MIDPOINT2)));
  static_assert(Warp::Midpoint4::INTERVALS ==
                (1U << static_cast<uint8_t>(CurlIntegrator::MIDPOINT4)));
  return static_cast<uint8_t>(1U << static_cast<uint8_t>(integrator));
}

inline constexpr const char *CURL_INTEGRATOR_IDS[] = {"euler-1", "midpoint-2",
                                                      "midpoint-4"};
static_assert(std::size(CURL_INTEGRATOR_IDS) ==
              static_cast<size_t>(CurlIntegrator::MIDPOINT4) + 1);

/** @brief Parameter family of warp.curl-flow.v2. */
struct CurlFlowWarpParams : Warp::CurlFlowParams {
  uint8_t basis = static_cast<uint8_t>(math::NoiseBasis::SIMPLEX);
  uint8_t integrator = static_cast<uint8_t>(CurlIntegrator::EULER1);
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<CurlFlowWarpParams>{
          "basis", &CurlFlowWarpParams::basis, NOISE_BASIS_IDS,
          static_cast<uint8_t>(math::NoiseBasis::SIMPLEX)},
      TopologyField<CurlFlowWarpParams>{"integrator",
                                        &CurlFlowWarpParams::integrator,
                                        CURL_INTEGRATOR_IDS, 0},
  };
};
static_assert(field_ids_unique<CurlFlowWarpParams>());
static_assert(
    appended_block_size_matches<CurlFlowWarpParams, Warp::CurlFlowParams, 2>(),
    "appended parameter block must have the expected rounded size");
static_assert(field_defaults_in_range<CurlFlowWarpParams>());

/** @brief The curl flow's prepared block: the owned noise field, this frame's
    point on the loop, and the sub-step count decoded from the integrator. */
struct PreparedCurlFlow {
  const FastNoiseLite *noise;
  math::Vector loop_offset;
  uint8_t intervals;
};

/** @brief PLANE endomorphism: flow along the component-clamped curl field. */
struct WarpCurlFlow : PhaseClockModel<NoisePhaseState> {
  static_assert([] {
    float max_scale = 0.0f;
    float max_strength = 0.0f;
    for (const auto &field : Warp::CurlFlowParams::FIELDS) {
      if (field.member == &Warp::CurlFlowParams::scale)
        max_scale = field.max;
      if (field.member == &Warp::CurlFlowParams::strength)
        max_strength = field.max > -field.min ? field.max : -field.min;
    }
    return max_scale * max_strength * Warp::CURL_VECTOR_COMPONENT_MAX <= 0.5f;
  }());
  static constexpr const char *ID = "warp.curl-flow.v2";
  static constexpr const char *NAME = "Curl Flow";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = CurlFlowWarpParams;
  using Prepared = PreparedCurlFlow;

  static void init(State &state, InstanceId id) { init_noise_phase(state, id); }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_noise_basis(params.basis);
    HS_CHECK(params.integrator <=
                 static_cast<uint8_t>(CurlIntegrator::MIDPOINT4),
             "warp.curl-flow: invalid integrator");
    const uint8_t intervals =
        curl_intervals(static_cast<CurlIntegrator>(params.integrator));
    return {&state.noise, math::noise_projected_loop_offset(state.phase),
            intervals};
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    if (!noise_plane_in_domain(input.coords))
      return input;
    return Kernel::warp(
        input, Warp::curl_flow(input.coords, *prepared.noise,
                               static_cast<math::NoiseBasis>(params.basis),
                               prepared.intervals, params.scale,
                               params.strength, prepared.loop_offset, true));
  }
};

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
