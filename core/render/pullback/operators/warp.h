/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "core/render/pullback/operators/common.h"
#include "core/render/pullback/authored_warp.h"

/**
 * @file warp.h
 * @brief PLANE-endomorphism operator models over the shared warp kernels.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Amplitude envelope of the warp operators. */
using WarpEnvelope = Warp::Envelope;

/** @brief Envelope topology values, in WarpEnvelope order. */
inline constexpr auto &WARP_ENVELOPE_IDS = Warp::ENVELOPE_IDS;

/**
 * @brief Bounds a warp operator's envelope enum8.
 * @details warp_envelope() relies on this check and carries no per-pixel guard.
 * @param envelope Raw `WarpEnvelope` value.
 */
inline void check_warp_envelope(uint8_t envelope) {
  HS_CHECK(envelope <= static_cast<uint8_t>(WarpEnvelope::EDGE_FADE),
           "warp operator: invalid envelope");
}

/** @brief The envelope switch over the shared warp envelope kernel.
 *  @param envelope Checked `WarpEnvelope` value.
 *  @param provenance Projection provenance of the plane sample.
 *  @param edge_width Edge-fade width.
 *  @return Amplitude factor for the sample. */
inline float warp_envelope(uint8_t envelope,
                           const ProjectionProvenance &provenance,
                           float edge_width) {
  return Warp::envelope(provenance, edge_width,
                        static_cast<WarpEnvelope>(envelope));
}

/** @brief Phase clock of the single-clock warp operators.
 * @details Column-time samples may lie outside [0, 1). */
struct WarpPhaseState {
  float phase = 0.0f; ///< Loop phase modulo 1.
};

/** @brief Prepares the clocks between consecutive endpoints.
 * @param from Previous state.
 * @param to Current state.
 * @return Prepared interval. */
inline Spacetime::ClockInterval temporal_interval(const WarpPhaseState &from,
                                                  const WarpPhaseState &to) {
  return {from.phase, to.phase, 1.0f};
}

/** @brief Applies the prepared clocks to a column state.
 * @tparam Params Operator parameter block.
 * @param out Copy of the current state.
 * @param interval Prepared frame interval.
 * @param fraction Frame fraction in [0, 1]. */
template <typename Params>
void temporal_sample(WarpPhaseState &out,
                     const Spacetime::ClockInterval &interval, const Params &,
                     float fraction) {
  out.phase = interval.sample(fraction);
}

/** @brief Path-length topology values: off, on. */
inline constexpr const char *AFFINE_PATH_LENGTH_IDS[] = {"off", "on"};

/** @brief Parameter family of warp.affine.v4: the affine controls plus whether
    the step contributes path length. */
struct AffineWarpParams : Warp::AffineParams {
  /** 1 measures the step's path length; 0 reports none. */
  uint8_t path_length = 1;

  static constexpr auto TOPOLOGY = std::array{
      TopologyField<AffineWarpParams>{"path-length",
                                      &AffineWarpParams::path_length,
                                      AFFINE_PATH_LENGTH_IDS, 1},
  };
};
static_assert(
    appended_block_size_matches<AffineWarpParams, Warp::AffineParams, 1>(),
    "appended parameter block must have the expected rounded size");

/** @brief Phase clock plus the accumulated frame rotation of warp.affine.v4.
 * @details Values are equivalent modulo their periods; column-time samples
 * may lie outside the canonical ranges. */
struct AffineClockState {
  float phase = 0.0f;    ///< Oscillation phase modulo 1.
  float rotation = 0.0f; ///< Frame rotation, radians modulo 2pi.
};

/** @brief Prepared frame clocks of AffineClockState. */
struct AffineClockInterval {
  Spacetime::ClockInterval phase;    ///< Prepared phase interval.
  Spacetime::ClockInterval rotation; ///< Prepared rotation interval.
};

/** @brief Prepares the clocks between consecutive endpoints.
 * @param from Previous state.
 * @param to Current state.
 * @return Prepared interval. */
inline AffineClockInterval temporal_interval(const AffineClockState &from,
                                             const AffineClockState &to) {
  return {{from.phase, to.phase, 1.0f},
          {from.rotation, to.rotation, math::TWO_PI_F}};
}

/** @brief Applies the prepared clocks to a column state.
 * @tparam Params Operator parameter block.
 * @param out Copy of the current state.
 * @param interval Prepared frame interval.
 * @param fraction Frame fraction in [0, 1]. */
template <typename Params>
void temporal_sample(AffineClockState &out, const AffineClockInterval &interval,
                     const Params &, float fraction) {
  out.phase = interval.phase.sample(fraction);
  out.rotation = interval.rotation.sample(fraction);
}

/** @brief PLANE endomorphism: the oscillating affine frame change, winding
    along the next stage's translation lattice. */
struct WarpAffineV4 : ValueStateModel<AffineClockState> {
  static constexpr const char *ID = "warp.affine.v4";
  static constexpr const char *NAME = "Affine Warp";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = AffineWarpParams;
  using Prepared = Warp::PreparedAffineSlot;

  static void advance(State &state, const Params &params) {
    state.phase = math::wrap_t(state.phase + params.speed);
    Warp::advance_affine_rotation(state.rotation, params);
  }
  /** @brief Bounds the warped coordinates over the affine cycle.
   * @param params Affine controls.
   * @param input_bound Maximum input plane radius.
   * @param lattice Downstream translation symmetry.
   * @return Maximum warped plane radius. */
  static float plane_bound(const Params &params, float input_bound,
                           const TranslationLattice &lattice) {
    return Warp::plane_bound(params, input_bound, lattice);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state,
                          const TranslationLattice &lattice) {
    return Warp::prepare(params, state.phase, state.rotation, lattice);
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    return Kernel::warp(input, Warp::affine_frame(input.coords, prepared,
                                                  params.path_length != 0));
  }
};

/** @brief Parameter family of warp.wave-shear.v2. */
using WaveShearWarpParams = Warp::AuthoredWaveShearParams;
static_assert(appended_block_size_matches<WaveShearWarpParams,
                                          Warp::WaveShearParams, 1>(),
              "appended parameter block must have the expected rounded size");

/** @brief The wave shear's prepared block: the field frame plus the phase the
    kernel consumes per sample. */
struct PreparedWaveShear {
  Warp::PreparedRotation rotation; ///< Shear field frame.
  float phase;                     ///< Normalized loop phase.
};

/** @brief PLANE endomorphism: the travelling sine shear. */
struct WarpWaveShear : PhaseClockModel<WarpPhaseState> {
  static constexpr const char *ID = "warp.wave-shear.v2";
  static constexpr const char *NAME = "Wave Shear";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = WaveShearWarpParams;
  using Prepared = PreparedWaveShear;

  /** @brief The factor covers the approximate sine.
   *  @param params Warp parameters.
   *  @param input_bound Bound on the input plane coordinates.
   *  @return Bound on the warped coordinates. */
  static float plane_bound(const Params &params, float input_bound) {
    return input_bound + 2.0f * fabsf(params.strength);
  }

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
    return Kernel::warp(
        input, Warp::wave_shear(input.coords, params, prepared.phase, amplitude,
                                prepared.rotation, MEASURE_PATH_LENGTH));
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

  /** @brief A rotation about a center no farther than the orbit's reach.
   *  @param params Warp parameters.
   *  @param input_bound Bound on the input plane coordinates.
   *  @return Bound on the warped coordinates. */
  static float plane_bound(const Params &params, float input_bound) {
    return input_bound + 2.0f * (std::hypot(params.center_x, params.center_y) +
                                 params.center_orbit_radius);
  }

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    return Warp::prepare(params, state.phase);
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &, const Prepared &prepared) {
    return Kernel::warp(
        input, Warp::vortex(input.coords, prepared, MEASURE_PATH_LENGTH));
  }
};

/** @brief Parameter family of warp.vector-noise.v2. */
using VectorNoiseWarpParams = Warp::AuthoredVectorNoiseParams;
static_assert(appended_block_size_matches<VectorNoiseWarpParams,
                                          Warp::VectorNoiseParams, 2>(),
              "appended parameter block must have the expected rounded size");

/** @brief The vector-noise warp's prepared block: the owned noise field plus
    the vector frame and loop point. */
struct PreparedVectorNoiseWarp {
  const FastNoiseLite *noise;         ///< The instance-owned noise field.
  Warp::PreparedVectorNoiseSlot slot; ///< Vector frame and loop point.
};

/** @brief PLANE endomorphism: the noise-vector displacement. */
struct WarpVectorNoise : PhaseClockModel<NoisePhaseState> {
  static constexpr const char *ID = "warp.vector-noise.v2";
  static constexpr const char *NAME = "Vector Noise";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = VectorNoiseWarpParams;
  using Prepared = PreparedVectorNoiseWarp;

  /** @brief The factor covers both noise channels with headroom.
   *  @param params Warp parameters.
   *  @param input_bound Bound on the input plane coordinates.
   *  @return Bound on the warped coordinates. */
  static float plane_bound(const Params &params, float input_bound) {
    return input_bound + 4.0f * fabsf(params.strength);
  }

  static void init(State &state, InstanceId id) { init_noise_phase(state, id); }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_warp_envelope(params.envelope);
    check_noise_basis(params.basis);
    return {&state.noise, Warp::prepare(params, state.phase)};
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    if (!math::noise_plane_in_domain(input.coords))
      return input;
    const float amplitude =
        params.strength *
        warp_envelope(params.envelope, input.provenance, params.edge_width);
    return Kernel::warp(
        input,
        Warp::vector_noise(input.coords, params, amplitude, *prepared.noise,
                           static_cast<math::NoiseBasis>(params.basis),
                           prepared.slot, MEASURE_PATH_LENGTH));
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

  /** @brief The fold lands inside one cell, whatever the input.
   *  @param params Warp parameters.
   *  @return The cell diagonal. */
  static float plane_bound(const Params &params, float) {
    return std::hypot(params.cell_x, params.cell_y);
  }

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    return Warp::prepare(params, state.phase);
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    return Kernel::warp(input, Warp::mirror_tile(input.coords, params, prepared,
                                                 MEASURE_PATH_LENGTH));
  }
};

using Warp::PolarMode;

/** @brief Parameter family of warp.polar-chart.v2. */
using PolarChartParams = Warp::AuthoredPolarParams;
static_assert(
    appended_block_size_matches<PolarChartParams, Warp::PolarParams, 2>(),
    "appended parameter block must have the expected rounded size");

/** @brief PLANE endomorphism: the polar chart change. */
struct WarpPolarChart : PhaseClockModel<WarpPhaseState> {
  static constexpr const char *ID = "warp.polar-chart.v2";
  static constexpr const char *NAME = "Polar Chart";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = PolarChartParams;
  /** @brief The frame's loop phase. */
  struct Prepared {
    float phase; ///< Normalized loop phase.
  };

  /** @brief Radial bound plus the angular reach of the harmonic, phases and
      atan2.
      @param params Warp parameters.
      @param input_bound Bound on the input plane coordinates.
      @return Bound on the warped coordinates. */
  static float plane_bound(const Params &params, float input_bound) {
    const float radial =
        static_cast<PolarMode>(params.mode) == PolarMode::LOGARITHMIC
            ? fmaxf(logf(fmaxf(input_bound, 1.0f)), logf(4096.0f))
            : input_bound;
    return params.radial_scale * radial + fabsf(params.radial_phase) +
           math::PI_F * static_cast<float>(params.harmonic + 2) +
           fabsf(params.angular_phase) + math::TWO_PI_F;
  }

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

using Warp::CurlIntegrator;
using Warp::curl_intervals;

/** @brief Parameter family of warp.curl-flow.v2. */
using CurlFlowWarpParams = Warp::AuthoredCurlFlowParams;
static_assert(
    appended_block_size_matches<CurlFlowWarpParams, Warp::CurlFlowParams, 2>(),
    "appended parameter block must have the expected rounded size");

/** @brief The curl flow's prepared block: the owned noise field, this frame's
    point on the loop, and the sub-step count decoded from the integrator. */
struct PreparedCurlFlow {
  const FastNoiseLite *noise; ///< The instance-owned noise field.
  math::Vector loop_offset;   ///< Lattice offset of this frame's loop phase.
  uint8_t intervals;          ///< Integrator sub-steps, curl_intervals().
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

  /** @brief Every sub-step follows a component-clamped vector.
   *  @param params Warp parameters.
   *  @param input_bound Bound on the input plane coordinates.
   *  @return Bound on the warped coordinates. */
  static float plane_bound(const Params &params, float input_bound) {
    return input_bound +
           fabsf(params.strength) * 2.0f * Warp::CURL_VECTOR_COMPONENT_MAX;
  }

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
    if (!math::noise_plane_in_domain(input.coords))
      return input;
    return Kernel::warp(
        input,
        Warp::curl_flow(input.coords, *prepared.noise,
                        static_cast<math::NoiseBasis>(params.basis),
                        prepared.intervals, params.scale, params.strength,
                        prepared.loop_offset, MEASURE_PATH_LENGTH));
  }
};

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
