/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file authored_warp.h
 * @brief Planar warp policies retaining authored topology controls. */

#include "core/render/pullback/warp.h"

namespace Pullback::Warp {

/** @brief Wave-shear controls with a runtime envelope. */
struct AuthoredWaveShearParams : WaveShearParams {
  uint8_t envelope = 0; ///< Envelope catalog index.
  /// Numeric controls, including runtime envelope width.
  static constexpr auto FIELDS = [] {
    auto fields = concat_fields<AuthoredWaveShearParams>(
        WaveShearParams::FIELDS,
        std::array<Field<AuthoredWaveShearParams>, 0>{});
    for (auto &field : fields)
      if (std::string_view(field.id) == "edge-width")
        field.name = "Warp Edge Width";
    return fields;
  }();
  /// Catalog values of the envelope selector.
  static constexpr auto &ENVELOPE_IDS = Warp::ENVELOPE_IDS;
  /// Authored topology selectors.
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<AuthoredWaveShearParams>{
          "envelope", &AuthoredWaveShearParams::envelope, ENVELOPE_IDS, 0},
  };
};
/** @brief Vector-noise controls with runtime basis and envelope. */
struct AuthoredVectorNoiseParams : VectorNoiseParams {
  uint8_t basis = 0;    ///< Noise-basis catalog index.
  uint8_t envelope = 0; ///< Envelope catalog index.
  /// Numeric controls, including runtime envelope width.
  static constexpr auto FIELDS = [] {
    auto fields = concat_fields<AuthoredVectorNoiseParams>(
        VectorNoiseParams::FIELDS,
        std::array<Field<AuthoredVectorNoiseParams>, 0>{});
    for (auto &field : fields)
      if (std::string_view(field.id) == "edge-width")
        field.name = "Warp Edge Width";
    return fields;
  }();
  /// Catalog values of the basis selector.
  static constexpr auto &BASIS_IDS = math::NOISE_BASIS_IDS;
  /// Catalog values of the envelope selector.
  static constexpr auto &ENVELOPE_IDS = Warp::ENVELOPE_IDS;
  /// Authored topology selectors.
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<AuthoredVectorNoiseParams>{
          "basis", &AuthoredVectorNoiseParams::basis, BASIS_IDS, 0},
      TopologyField<AuthoredVectorNoiseParams>{
          "envelope", &AuthoredVectorNoiseParams::envelope, ENVELOPE_IDS, 0},
  };
};
/** @brief Polar chart controls with runtime mode and harmonic. */
struct AuthoredPolarParams : PolarParams {
  uint8_t mode = 0;     ///< Linear or logarithmic catalog index.
  uint8_t harmonic = 0; ///< Harmonic minus one.
  /// Catalog values of the mode selector.
  static constexpr const char *MODE_IDS[] = {"linear", "logarithmic"};
  /// Catalog values of the harmonic selector.
  static constexpr const char *HARMONIC_IDS[] = {
      "h1", "h2",  "h3",  "h4",  "h5",  "h6",  "h7",  "h8",
      "h9", "h10", "h11", "h12", "h13", "h14", "h15", "h16"};
  /// Authored topology selectors.
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<AuthoredPolarParams>{"mode", &AuthoredPolarParams::mode,
                                         MODE_IDS, 0},
      TopologyField<AuthoredPolarParams>{
          "harmonic", &AuthoredPolarParams::harmonic, HARMONIC_IDS, 0},
  };
};
/** @brief Curl-flow controls with runtime basis and integrator. */
struct AuthoredCurlFlowParams : CurlFlowParams {
  uint8_t basis = 0;      ///< Noise-basis catalog index.
  uint8_t integrator = 0; ///< Integration substeps are two to this power.
  /// Catalog values of the basis selector.
  static constexpr auto &BASIS_IDS = math::NOISE_BASIS_IDS;
  /// Catalog values of the integrator selector.
  static constexpr const char *INTEGRATOR_IDS[] = {"euler-1", "midpoint-2",
                                                   "midpoint-4"};
  /// Authored topology selectors.
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<AuthoredCurlFlowParams>{
          "basis", &AuthoredCurlFlowParams::basis, BASIS_IDS, 0},
      TopologyField<AuthoredCurlFlowParams>{
          "integrator", &AuthoredCurlFlowParams::integrator, INTEGRATOR_IDS, 0},
  };
};

/** @brief Wave shear with a runtime-selected amplitude envelope.
 * @tparam State Keyed warp provider. */
template <typename State> struct ConfiguredWaveShear : WaveShear<State> {
  using Base = WaveShear<State>;                ///< Static policy contract.
  using FrameState = typename Base::FrameState; ///< Bound frame.
  using Prepared = typename Base::Prepared;     ///< Vector rotation.
  /** @brief Prepares and validates the envelope.
   * @param frame Current frame.
   * @return Vector rotation. */
  static Prepared prepare(const FrameState &frame) {
    HS_CHECK(State::params(frame).envelope <=
                 static_cast<uint8_t>(Envelope::EDGE_FADE),
             "invalid warp envelope");
    return State::prepare(frame);
  }
  /** @brief Applies the authored sine shear.
   * @param input Projected coordinate.
   * @param provenance Projection metadata.
   * @param frame Current frame.
   * @param prepared Vector rotation.
   * @return Warped coordinate and path length. */
  static WarpStepResult apply(const math::Complex &input,
                              const ProjectionProvenance &provenance,
                              const FrameState &frame,
                              const Prepared &prepared) {
    const auto &params = State::params(frame);
    return wave_shear(input, params, State::phase(frame),
                      params.strength *
                          envelope(provenance, params.edge_width,
                                   static_cast<Envelope>(params.envelope)),
                      prepared, State::path_length_required(frame));
  }
};

/** @brief Vector noise with runtime basis and envelope.
 * @tparam State Keyed warp provider. */
template <typename State>
struct ConfiguredVectorNoise
    : VectorNoise<State, math::NoiseBasis::SIMPLEX, FlatEnvelope> {
  using Base = VectorNoise<State, math::NoiseBasis::SIMPLEX,
                           FlatEnvelope>;       ///< Static policy contract.
  using FrameState = typename Base::FrameState; ///< Bound frame.
  using Prepared = typename Base::Prepared; ///< Noise-loop and vector frame.
  /** @brief Prepares and validates the authored topology.
   * @param frame Current frame.
   * @return Noise-loop and vector frame. */
  static Prepared prepare(const FrameState &frame) {
    const auto &params = State::params(frame);
    HS_CHECK(params.basis <= static_cast<uint8_t>(math::NoiseBasis::RIDGED3) &&
                 params.envelope <= static_cast<uint8_t>(Envelope::EDGE_FADE),
             "invalid vector-noise topology");
    return State::prepare(frame);
  }
  /** @brief Applies the authored noise-vector displacement.
   * @param input Projected coordinate.
   * @param provenance Projection metadata.
   * @param frame Current frame.
   * @param prepared Noise-loop and vector frame.
   * @return Warped coordinate and path length. */
  static WarpStepResult apply(const math::Complex &input,
                              const ProjectionProvenance &provenance,
                              const FrameState &frame,
                              const Prepared &prepared) {
    if (!math::noise_plane_in_domain(input))
      return {input, 0.0f};
    const auto &params = State::params(frame);
    return vector_noise(
        input, params,
        params.strength * envelope(provenance, params.edge_width,
                                   static_cast<Envelope>(params.envelope)),
        State::noise(frame), static_cast<math::NoiseBasis>(params.basis),
        prepared, State::path_length_required(frame));
  }
};

/** @brief Polar chart with runtime mode and harmonic.
 * @tparam State Keyed warp provider. */
template <typename State>
struct ConfiguredPolarChart : PolarChart<State, LinearPolar, 1> {
  using Prepared = NoPrepared;                   ///< No per-frame preparation.
  using FrameState = typename State::FrameState; ///< Bound frame.
  /** @brief Validates the mode and harmonic.
   * @param frame Current frame.
   * @return Empty preparation. */
  static NoPrepared prepare(const FrameState &frame) {
    const auto &params = State::params(frame);
    HS_CHECK(params.mode < std::size(AuthoredPolarParams::MODE_IDS) &&
                 params.harmonic < std::size(AuthoredPolarParams::HARMONIC_IDS),
             "invalid polar-chart topology");
    return {};
  }
  /** @brief Applies the authored polar chart.
   * @param input Projected coordinate.
   * @param frame Current frame.
   * @return Warped coordinate and zero path length. */
  static WarpStepResult apply(const math::Complex &input,
                              const ProjectionProvenance &,
                              const FrameState &frame, const NoPrepared &) {
    const auto &params = State::params(frame);
    return polar_chart(input, params, State::phase(frame), params.mode == 1,
                       static_cast<uint8_t>(params.harmonic + 1));
  }
};

/** @brief Curl flow with runtime basis and integrator.
 * @tparam State Keyed warp provider. */
template <typename State>
struct ConfiguredCurlFlow : CurlFlow<State, math::NoiseBasis::SIMPLEX, Euler1> {
  using Base = CurlFlow<State, math::NoiseBasis::SIMPLEX,
                        Euler1>;                ///< Static policy contract.
  using FrameState = typename Base::FrameState; ///< Bound frame.
  using Prepared = typename Base::Prepared;     ///< Noise-loop offset.
  /** @brief Prepares and validates the authored topology.
   * @param frame Current frame.
   * @return Noise-loop offset. */
  static Prepared prepare(const FrameState &frame) {
    const auto &params = State::params(frame);
    HS_CHECK(params.basis <= static_cast<uint8_t>(math::NoiseBasis::RIDGED3) &&
                 params.integrator <
                     std::size(AuthoredCurlFlowParams::INTEGRATOR_IDS),
             "invalid curl-flow topology");
    return Base::prepare(frame);
  }
  /** @brief Applies the authored curl flow.
   * @param input Projected coordinate.
   * @param frame Current frame.
   * @param prepared Noise-loop offset.
   * @return Warped coordinate and path length. */
  static WarpStepResult apply(const math::Complex &input,
                              const ProjectionProvenance &,
                              const FrameState &frame,
                              const Prepared &prepared) {
    if (!math::noise_plane_in_domain(input))
      return {input, 0.0f};
    const auto &params = State::params(frame);
    return curl_flow(
        input, State::noise(frame), static_cast<math::NoiseBasis>(params.basis),
        static_cast<uint8_t>(1u << params.integrator), params.scale,
        params.strength, prepared, State::path_length_required(frame));
  }
};

} // namespace Pullback::Warp
