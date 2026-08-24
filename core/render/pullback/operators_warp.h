/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "render/pullback/operators_common.h"
#include "render/pullback/warp.h"

/**
 * @file operators_warp.h
 * @brief PLANE-endomorphism operator models over the shared warp kernels.
 */

namespace Pullback {

namespace Interp {

namespace Op {

enum class WarpEnvelope : uint8_t {
  FLAT = 0,
  PROJECTION_WEIGHT = 1,
  EDGE_FADE = 2
};

inline constexpr const char *WARP_ENVELOPE_IDS[] = {"flat", "projection-weight",
                                                    "edge-fade"};

/** @brief The envelope switch over the shared warp envelope kernel. */
inline float warp_envelope(uint8_t envelope,
                           const ProjectionProvenance &provenance,
                           float edge_width) {
  const auto mode = static_cast<WarpEnvelope>(envelope);
  return Warp::envelope(provenance, edge_width,
                        mode == WarpEnvelope::PROJECTION_WEIGHT,
                        mode == WarpEnvelope::EDGE_FADE);
}

/** @brief Phase clock of the single-clock warp operators. */
struct WarpPhaseState {
  float phase = 0.0f;
};

/** @brief Parameter family of warp.affine.v2.
    @details Translation is in plane units per phase turn: the chain has no
    lattice source coupling, so the composed path's cell scaling is fixed
    at 1. */
struct AffineWarpParams : Warp::AffineParams {
  static constexpr auto FIELDS = concat_fields<AffineWarpParams>(
      Warp::AffineParams::FIELDS, std::array<Field<AffineWarpParams>, 0>{});
};
static_assert(field_ids_unique<AffineWarpParams>());

/** @brief Phase clock plus the accumulated frame rotation of warp.affine.v2. */
struct AffineClockState {
  float phase = 0.0f;
  float rotation = 0.0f;
};

/** @brief PLANE endomorphism: the oscillating affine frame change. */
struct WarpAffine {
  static constexpr const char *ID = "warp.affine.v2";
  static constexpr const char *NAME = "Affine Warp";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = AffineWarpParams;
  using State = AffineClockState;
  using Prepared = Warp::PreparedAffineSlot;

  static void init(State &, InstanceId) {}
  static Status migrate(State &dst, const State &src, InstanceId) {
    dst = src;
    return Status::OK;
  }
  static void advance(State &state, const Params &params) {
    state.phase = wrap_t(state.phase + params.speed);
    state.rotation =
        TWO_PI_F *
        wrap_t((state.rotation + params.speed * params.rotation_rate) /
               TWO_PI_F);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    return Warp::prepare(params, state.phase, state.rotation, 1.0f);
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

  static constexpr auto FIELDS = concat_fields<WaveShearWarpParams>(
      Warp::WaveShearParams::FIELDS,
      std::array<Field<WaveShearWarpParams>, 0>{});
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<WaveShearWarpParams>{
          "envelope", &WaveShearWarpParams::envelope, WARP_ENVELOPE_IDS, 3,
          static_cast<uint8_t>(WarpEnvelope::FLAT)},
  };
};
static_assert(field_ids_unique<WaveShearWarpParams>());

/** @brief The wave shear's prepared block: the field frame plus the phase the
    kernel consumes per sample. */
struct PreparedWaveShear {
  Warp::PreparedRotation rotation;
  float phase;
};

/** @brief PLANE endomorphism: the travelling sine shear. */
struct WarpWaveShear {
  static constexpr const char *ID = "warp.wave-shear.v2";
  static constexpr const char *NAME = "Wave Shear";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = WaveShearWarpParams;
  using State = WarpPhaseState;
  using Prepared = PreparedWaveShear;

  static void init(State &, InstanceId) {}
  static Status migrate(State &dst, const State &src, InstanceId) {
    dst = src;
    return Status::OK;
  }
  static void advance(State &state, const Params &params) {
    state.phase = wrap_t(state.phase + params.speed);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
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
struct VortexWarpParams {
  float speed = 0.0f;
  float center_x = 0.0f;
  float center_y = 0.0f;
  float radius = 1.0f;
  float turns = 0.0f;
  float center_orbit_radius = 0.0f;

  static constexpr auto FIELDS = std::array{
      Field<VortexWarpParams>{"speed", &VortexWarpParams::speed, nullptr,
                              -0.02f, 0.02f, FieldCurve::LERP},
      Field<VortexWarpParams>{"center-x", &VortexWarpParams::center_x,
                              "Vortex Center X", -4.0f, 4.0f, FieldCurve::LERP},
      Field<VortexWarpParams>{"center-y", &VortexWarpParams::center_y,
                              "Vortex Center Y", -4.0f, 4.0f, FieldCurve::LERP},
      Field<VortexWarpParams>{"radius", &VortexWarpParams::radius,
                              "Vortex Radius", 1.0f / 64.0f, 8.0f,
                              FieldCurve::LOG_POSITIVE},
      Field<VortexWarpParams>{"turns", &VortexWarpParams::turns, "Vortex Turns",
                              -4.0f, 4.0f, FieldCurve::LERP},
      Field<VortexWarpParams>{
          "center-orbit-radius", &VortexWarpParams::center_orbit_radius,
          "Vortex Center Orbit", 0.0f, 4.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<VortexWarpParams>());

struct PreparedVortexWarp {
  struct {
    struct {
      float center_x;
      float center_y;
      float radius_sq;
      float angle_numerator;
    } vortex;
  } transform;
};

/** @brief PLANE endomorphism: the orbiting radial vortex. */
struct WarpVortex {
  static constexpr const char *ID = "warp.vortex.v2";
  static constexpr const char *NAME = "Vortex";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = VortexWarpParams;
  using State = WarpPhaseState;
  using Prepared = PreparedVortexWarp;

  static void init(State &, InstanceId) {}
  static Status migrate(State &dst, const State &src, InstanceId) {
    dst = src;
    return Status::OK;
  }
  static void advance(State &state, const Params &params) {
    state.phase = wrap_t(state.phase + params.speed);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    const float phase = TWO_PI_F * state.phase;
    return {{{params.center_x + params.center_orbit_radius * cosf(phase),
              params.center_y + params.center_orbit_radius * sinf(phase),
              params.radius * params.radius, TWO_PI_F * params.turns}}};
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &, const Prepared &prepared) {
    return Kernel::warp(input, Warp::vortex(input.coords, prepared, true));
  }
};

/** @brief Parameter family of warp.vector-noise.v2. */
struct VectorNoiseWarpParams : Warp::VectorNoiseParams {
  uint8_t basis = static_cast<uint8_t>(::NoiseBasis::SIMPLEX);
  uint8_t envelope = static_cast<uint8_t>(WarpEnvelope::FLAT);

  static constexpr auto FIELDS = concat_fields<VectorNoiseWarpParams>(
      Warp::VectorNoiseParams::FIELDS,
      std::array<Field<VectorNoiseWarpParams>, 0>{});
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<VectorNoiseWarpParams>{
          "basis", &VectorNoiseWarpParams::basis, NOISE_BASIS_IDS, 3,
          static_cast<uint8_t>(::NoiseBasis::SIMPLEX)},
      TopologyField<VectorNoiseWarpParams>{
          "envelope", &VectorNoiseWarpParams::envelope, WARP_ENVELOPE_IDS, 3,
          static_cast<uint8_t>(WarpEnvelope::FLAT)},
  };
};
static_assert(field_ids_unique<VectorNoiseWarpParams>());

/** @brief The vector-noise warp's prepared block: the owned noise field plus
    the vector frame and loop point. */
struct PreparedVectorNoiseWarp {
  const FastNoiseLite *noise;
  Warp::PreparedVectorNoiseSlot slot;
};

/** @brief PLANE endomorphism: the noise-vector displacement. */
struct WarpVectorNoise {
  static constexpr const char *ID = "warp.vector-noise.v2";
  static constexpr const char *NAME = "Vector Noise";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = VectorNoiseWarpParams;
  using State = NoisePhaseState;
  using Prepared = PreparedVectorNoiseWarp;

  static void init(State &state, InstanceId id) { init_noise_phase(state, id); }
  static Status migrate(State &dst, const State &src, InstanceId) {
    dst = src;
    return Status::OK;
  }
  static void advance(State &state, const Params &params) {
    state.phase = wrap_t(state.phase + params.speed);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    return {&state.noise, Warp::prepare(params, state.phase)};
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    const float amplitude =
        params.strength *
        warp_envelope(params.envelope, input.provenance, params.edge_width);
    return Kernel::warp(
        input,
        Warp::vector_noise(input.coords, params, amplitude, *prepared.noise,
                           static_cast<::NoiseBasis>(params.basis),
                           prepared.slot, true));
  }
};

/** @brief Parameter family of warp.mirror-tile.v2. */
struct MirrorWarpParams : Warp::MirrorParams {
  static constexpr auto FIELDS = concat_fields<MirrorWarpParams>(
      Warp::MirrorParams::FIELDS, std::array<Field<MirrorWarpParams>, 0>{});
};
static_assert(field_ids_unique<MirrorWarpParams>());

/** @brief PLANE endomorphism: the mirrored tiling fold. */
struct WarpMirrorTile {
  static constexpr const char *ID = "warp.mirror-tile.v2";
  static constexpr const char *NAME = "Mirror Tile";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = MirrorWarpParams;
  using State = WarpPhaseState;
  using Prepared = Warp::PreparedMirrorSlot;

  static void init(State &, InstanceId) {}
  static Status migrate(State &dst, const State &src, InstanceId) {
    dst = src;
    return Status::OK;
  }
  static void advance(State &state, const Params &params) {
    state.phase = wrap_t(state.phase + params.speed);
  }
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

/** @brief Harmonic topology values h1..h16, indexed by harmonic - 1. */
inline constexpr const char *POLAR_HARMONIC_IDS[] = {
    "h1", "h2",  "h3",  "h4",  "h5",  "h6",  "h7",  "h8",
    "h9", "h10", "h11", "h12", "h13", "h14", "h15", "h16"};
static_assert(std::size(POLAR_HARMONIC_IDS) == Warp::MAX_POLAR_HARMONIC);

/** @brief Parameter family of warp.polar-chart.v2. */
struct PolarChartParams : Warp::PolarParams {
  uint8_t mode = static_cast<uint8_t>(PolarMode::LINEAR);
  uint8_t harmonic = 0; /**< Harmonic value index; harmonic = index + 1. */

  static constexpr auto FIELDS = concat_fields<PolarChartParams>(
      Warp::PolarParams::FIELDS, std::array<Field<PolarChartParams>, 0>{});
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<PolarChartParams>{"mode", &PolarChartParams::mode,
                                      POLAR_MODE_IDS, 2,
                                      static_cast<uint8_t>(PolarMode::LINEAR)},
      TopologyField<PolarChartParams>{"harmonic", &PolarChartParams::harmonic,
                                      POLAR_HARMONIC_IDS,
                                      Warp::MAX_POLAR_HARMONIC, 0},
  };
};
static_assert(field_ids_unique<PolarChartParams>());

/** @brief PLANE endomorphism: the polar chart change. */
struct WarpPolarChart {
  static constexpr const char *ID = "warp.polar-chart.v2";
  static constexpr const char *NAME = "Polar Chart";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = PolarChartParams;
  using State = WarpPhaseState;
  struct Prepared {
    float phase;
  };

  static void init(State &, InstanceId) {}
  static Status migrate(State &dst, const State &src, InstanceId) {
    dst = src;
    return Status::OK;
  }
  static void advance(State &state, const Params &params) {
    state.phase = wrap_t(state.phase + params.speed);
  }
  static Prepared prepare(const FrameContext &, const Params &,
                          const State &state) {
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

inline constexpr const char *CURL_INTEGRATOR_IDS[] = {"euler-1", "midpoint-2",
                                                      "midpoint-4"};

/** @brief Parameter family of warp.curl-flow.v2. */
struct CurlFlowParams {
  float speed = 0.0f;    /**< Per-frame advance of the flow's loop phase. */
  float strength = 0.0f; /**< Flow distance; 0 skips the stage. */
  float scale = 1.0f;    /**< Spatial scale of the sampled field. */
  uint8_t basis = static_cast<uint8_t>(::NoiseBasis::SIMPLEX);
  uint8_t integrator = 0;

  static constexpr auto FIELDS = std::array{
      Field<CurlFlowParams>{"speed", &CurlFlowParams::speed, nullptr, -0.02f,
                            0.02f, FieldCurve::LERP},
      Field<CurlFlowParams>{"strength", &CurlFlowParams::strength,
                            "Warp Strength", -30.0f, 30.0f, FieldCurve::LERP},
      Field<CurlFlowParams>{"scale", &CurlFlowParams::scale, "Warp Scale",
                            1.0f / 64.0f, 64.0f, FieldCurve::LOG_POSITIVE},
  };
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<CurlFlowParams>{
          "basis", &CurlFlowParams::basis, NOISE_BASIS_IDS, 3,
          static_cast<uint8_t>(::NoiseBasis::SIMPLEX)},
      TopologyField<CurlFlowParams>{"integrator", &CurlFlowParams::integrator,
                                    CURL_INTEGRATOR_IDS, 3, 0},
  };
};
static_assert(field_ids_unique<CurlFlowParams>());

/** @brief The curl flow's prepared block: the owned noise field plus the loop
    phase. */
struct PreparedCurlFlow {
  const FastNoiseLite *noise;
  float phase;
  uint8_t intervals;
};

/** @brief PLANE endomorphism: the divergence-free curl flow, on the
    single-interval integrator. */
struct WarpCurlFlow {
  static constexpr const char *ID = "warp.curl-flow.v2";
  static constexpr const char *NAME = "Curl Flow";
  using Input = PlaneSample;
  using Output = PlaneSample;
  using Params = CurlFlowParams;
  using State = NoisePhaseState;
  using Prepared = PreparedCurlFlow;

  static void init(State &state, InstanceId id) { init_noise_phase(state, id); }
  static Status migrate(State &dst, const State &src, InstanceId) {
    dst = src;
    return Status::OK;
  }
  static void advance(State &state, const Params &params) {
    state.phase = wrap_t(state.phase + params.speed);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    return {&state.noise, state.phase,
            static_cast<uint8_t>(1U << params.integrator)};
  }
  static PlaneSample run(const PlaneSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    return Kernel::warp(input,
                        Warp::curl_flow(input.coords, *prepared.noise,
                                        static_cast<::NoiseBasis>(params.basis),
                                        prepared.intervals, params.scale,
                                        params.strength, prepared.phase, true));
  }
};

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
