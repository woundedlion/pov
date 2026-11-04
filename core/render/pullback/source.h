/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>

#include "core/pattern/noise.h"
#include "core/pattern/sources.h"

#include "core/math/projection_patterns.h"
#include "core/render/pullback/contract.h"
#include "core/render/pullback/fields.h"
#include "core/render/pullback/translation_lattice.h"
#include "core/math/3dmath.h"

/**
 * @file source.h
 * @brief Scalar source-field policies.
 */

namespace Pullback {

namespace Source {

/** @brief Borrows canonical controls for a shared source kernel.
 * @tparam Base Canonical source controls.
 * @tparam Params Provider controls.
 * @param params Provider parameter family.
 * @return Canonical base when available, otherwise the provider controls. */
template <typename Base, typename Params>
HS_HOT_INLINE const auto &source_kernel_params(const Params &params) {
  if constexpr (std::derived_from<Params, Base>)
    return static_cast<const Base &>(params);
  else
    return params;
}

/** @brief Parameters of the concentric ring source. */
struct RingsSourceParams : Pattern::RingsSourceParams {
  /** @brief Per-parameter id, member, label, range and curve. */
  static constexpr auto FIELDS = std::array{
      Field<RingsSourceParams>{"pattern-freq", &RingsSourceParams::pattern_freq,
                               "Pattern Freq", 0.1f, 20.0f, FieldCurve::LERP},
      Field<RingsSourceParams>{"speed", &RingsSourceParams::speed, "Speed",
                               0.0f, 0.5f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<RingsSourceParams>());
static_assert(field_defaults_in_range<RingsSourceParams>());

/**
 * @brief Source parameters for the coupled sine grid
 *        (Pullback::Source::Grid).
 */
struct GridSourceParams : Pattern::GridSourceParams {
  /** @brief Per-parameter id, member, label, range and curve. */
  static constexpr auto FIELDS = std::array{
      Field<GridSourceParams>{"pattern-freq", &GridSourceParams::pattern_freq,
                              "Pattern Freq", 0.01f, 64.0f, FieldCurve::LERP},
      Field<GridSourceParams>{"speed", &GridSourceParams::speed, "Speed", 0.0f,
                              0.5f, FieldCurve::LERP},
      Field<GridSourceParams>{"complexity", &GridSourceParams::complexity,
                              "Complexity", 0.0f, 3.0f, FieldCurve::LERP},
      Field<GridSourceParams>{"pattern-mix", &GridSourceParams::pattern_mix,
                              "Pattern Mix", 0.0f, 1.0f, FieldCurve::LERP},
      Field<GridSourceParams>{"drift", &GridSourceParams::secondary_rate,
                              "Drift", 0.0f, 1.25f, FieldCurve::LERP},
      Field<GridSourceParams>{"angle-speed", &GridSourceParams::angle_rate,
                              "Source Angle Speed", 0.0f, 0.05f,
                              FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<GridSourceParams>());
static_assert(field_defaults_in_range<GridSourceParams>());

/**
 * @brief Source parameters for the two-wave interference field
 *        (Pullback::Source::TwinWave).
 */
struct TwinWaveSourceParams : Pattern::TwinWaveSourceParams {
  /** @brief Per-parameter id, member, label, range and curve. */
  static constexpr auto FIELDS = std::array{
      Field<TwinWaveSourceParams>{
          "pattern-freq", &TwinWaveSourceParams::pattern_freq, "Pattern Freq",
          0.1f, 20.0f, FieldCurve::LERP},
      Field<TwinWaveSourceParams>{"speed", &TwinWaveSourceParams::speed,
                                  "Speed", 0.0f, 0.5f, FieldCurve::LERP},
      Field<TwinWaveSourceParams>{"drift",
                                  &TwinWaveSourceParams::secondary_rate,
                                  "Drift", 0.0f, 1.25f, FieldCurve::LERP},
      Field<TwinWaveSourceParams>{
          "angle-speed", &TwinWaveSourceParams::angle_rate,
          "Source Angle Speed", 0.0f, 0.05f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<TwinWaveSourceParams>());
static_assert(field_defaults_in_range<TwinWaveSourceParams>());

/**
 * @brief Source parameters for the rotating spiral field
 *        (Pullback::Source::Spiral).
 */
struct SpiralSourceParams : Pattern::SpiralSourceParams {
  /** @brief Per-parameter id, member, label, range and curve. */
  static constexpr auto FIELDS = std::array{
      Field<SpiralSourceParams>{"pattern-freq",
                                &SpiralSourceParams::pattern_freq,
                                "Pattern Freq", 0.1f, 20.0f, FieldCurve::LERP},
      Field<SpiralSourceParams>{"speed", &SpiralSourceParams::speed, "Speed",
                                0.0f, 0.5f, FieldCurve::LERP},
      Field<SpiralSourceParams>{"angle-speed", &SpiralSourceParams::angle_rate,
                                "Source Angle Speed", 0.0f, 0.05f,
                                FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<SpiralSourceParams>());
static_assert(field_defaults_in_range<SpiralSourceParams>());

/**
 * @brief Source parameters for the noise-contour sources
 *        (Pullback::Source::ProjectedNoise).
 */
struct NoiseSourceParams : Pattern::NoiseSourceParams {
  /** @brief Per-parameter id, member, label, range and curve. */
  static constexpr auto FIELDS = std::array{
      Field<NoiseSourceParams>{"noise-scale", &NoiseSourceParams::noise_scale,
                               "Source Noise Scale", 1.0f / 64.0f, 64.0f,
                               FieldCurve::LOG_POSITIVE},
      Field<NoiseSourceParams>{
          "noise-contrast", &NoiseSourceParams::noise_contrast,
          "Source Noise Contrast", 0.0f, 8.0f, FieldCurve::LERP},
      Field<NoiseSourceParams>{
          "noise-speed", &NoiseSourceParams::noise_time_rate,
          "Source Noise Speed", -1.0f / 64.0f, 1.0f / 64.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<NoiseSourceParams>());
static_assert(field_defaults_in_range<NoiseSourceParams>());

/**
 * @brief Source family selecting the plane-domain noise contour
 *        (Pullback::Source::ProjectedNoise): the field is sampled in the
 *        projected chart, so it carries the projection's distortion.
 */
struct ProjectedNoiseSourceParams : NoiseSourceParams {};

/**
 * @brief Source parameters for the per-cell primitive lattice
 *        (Pullback::Source::PrimitiveLattice).
 */
struct LatticeSourceParams : Pattern::LatticeSourceParams {
  /// translation_lattice() is axis-aligned for every parameter value.
  static constexpr bool RECTANGULAR_TRANSLATION_LATTICE = true;
  /** @brief Per-parameter id, member, label, range and curve. */
  static constexpr auto FIELDS = std::array{
      Field<LatticeSourceParams>{
          "lattice-cell-scale", &LatticeSourceParams::lattice_cell_scale,
          "Lattice Cell Scale", 1.0f / 64.0f, 8.0f, FieldCurve::LOG_POSITIVE},
      Field<LatticeSourceParams>{"lattice-shape",
                                 &LatticeSourceParams::lattice_shape_blend,
                                 "Lattice Shape", 0.0f, 1.0f, FieldCurve::LERP},
      Field<LatticeSourceParams>{
          "lattice-softness", &LatticeSourceParams::lattice_softness,
          "Lattice Softness", 1.0f / 1024.0f, 1.0f, FieldCurve::LOG_POSITIVE},
      Field<LatticeSourceParams>{
          "lattice-radius", &LatticeSourceParams::lattice_radius,
          "Lattice Radius", 1.0f / 64.0f, 0.49f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<LatticeSourceParams>());
static_assert(field_defaults_in_range<LatticeSourceParams>());

/** @brief Square translation symmetry of the primitive lattice.
    @param params Lattice controls.
    @return Square lattice of period `1 / lattice_cell_scale`. */
constexpr TranslationLattice
translation_lattice(const LatticeSourceParams &params) {
  return square_lattice(1.0f / params.lattice_cell_scale);
}

/** @brief Source parameters for latitude bands on a moving sphere. */
/** @brief Source parameters for periodic polygon edge tessellations. */
struct TessellationSourceParams : Pattern::TessellationSourceParams {
  /** @brief Per-parameter id, member, label, range and curve. */
  static constexpr auto FIELDS = std::array{
      Field<TessellationSourceParams>{
          "cell-scale", &TessellationSourceParams::cell_scale, "Cell Scale",
          1.0f / 64.0f, 8.0f, FieldCurve::LOG_POSITIVE},
      Field<TessellationSourceParams>{
          "line-thickness", &TessellationSourceParams::line_thickness,
          "Line Thickness", 1.0f / 1024.0f, 0.25f, FieldCurve::LOG_POSITIVE},
      Field<TessellationSourceParams>{
          "line-softness", &TessellationSourceParams::line_softness,
          "Line Softness", 1.0f / 1024.0f, 0.25f, FieldCurve::LOG_POSITIVE},
      Field<TessellationSourceParams>{
          "angle-speed", &TessellationSourceParams::angle_rate,
          "Tessellation Spin Speed", -0.05f, 0.05f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<TessellationSourceParams>());
static_assert(field_defaults_in_range<TessellationSourceParams>());

/** @brief A parameter provider whose `prepare` builds the source's per-frame
    state. */
template <typename State, typename Binding>
concept StateProvider = Detail::ParamsProvider<State, Binding> &&
                        requires(const typename Binding::FrameState &frame) {
                          State::prepare(frame);
                        };

/**
 * @brief TwinWave source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), prepare(frame) accessors.
 */
template <typename State> struct TwinWave : PrepareFromState<State> {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the pattern
   *        frequency, rotation and both phases.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::prepare(frame).angle_cos } -> std::convertible_to<float>;
        { State::prepare(frame).angle_sin } -> std::convertible_to<float>;
        { State::prepare(frame).primary } -> std::convertible_to<float>;
        { State::prepare(frame).secondary } -> std::convertible_to<float>;
      };

  /// Phases `State::prepare` resolves for the frame.
  using Prepared = StatePrepared<State>;

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return Pattern::twin_wave(
        Pattern::PlanePosition{projections::stereo_pattern_args(
            input.coords, params.pattern_freq)},
        prepared);
  }
};

/**
 * @brief Rings source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), prepare(frame) accessors.
 */
template <typename State> struct Rings : PrepareFromState<State> {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the pattern
   *        frequency and primary phase.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::prepare(frame).primary } -> std::convertible_to<float>;
      };

  /// Phases `State::prepare` resolves for the frame.
  using Prepared = StatePrepared<State>;

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return Pattern::rings(
        Pattern::PlanePosition{projections::stereo_pattern_args(
            input.coords, params.pattern_freq)},
        prepared);
  }
};

/**
 * @brief Spiral source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), prepare(frame) accessors.
 */
template <typename State> struct Spiral : PrepareFromState<State> {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the pattern
   *        frequency, angle and primary phase.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::prepare(frame).angle } -> std::convertible_to<float>;
        { State::prepare(frame).primary } -> std::convertible_to<float>;
      };

  /// Phases `State::prepare` resolves for the frame.
  using Prepared = StatePrepared<State>;

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return Pattern::spiral(
        Pattern::PlanePosition{projections::stereo_pattern_args(
            input.coords, params.pattern_freq)},
        prepared);
  }
};

/**
 * @brief Grid source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), prepare(frame) accessors.
 */
template <typename State> struct Grid : PrepareFromState<State> {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the grid
   *        parameters, rotation and both phases.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::params(frame).pattern_mix } -> std::convertible_to<float>;
        { State::params(frame).complexity } -> std::convertible_to<float>;
        { State::prepare(frame).angle_cos } -> std::convertible_to<float>;
        { State::prepare(frame).angle_sin } -> std::convertible_to<float>;
        { State::prepare(frame).primary } -> std::convertible_to<float>;
        { State::prepare(frame).secondary } -> std::convertible_to<float>;
      };

  /// Phase carrier returned by the provider.
  using Phases = StatePrepared<State>;
  /// Prepared phases with optional stationary kernel selection.
  using Prepared =
      std::conditional_t<std::is_same_v<Phases, Pattern::PreparedSource>,
                         Pattern::PreparedGrid, Phases>;

  /** @brief Selects the grid kernel from this column's phases and parameters.
   * @param frame Column frame state.
   * @return Prepared phases and kernel selection. */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    if constexpr (std::is_same_v<Prepared, Pattern::PreparedGrid>) {
      const auto phases = State::prepare(frame);
      const auto &params = State::params(frame);
      return {phases, params.complexity == 0.0f && params.pattern_mix == 0.0f &&
                          phases.primary == 0.0f && phases.secondary == 0.0f &&
                          phases.angle_cos == 1.0f && phases.angle_sin == 0.0f};
    } else {
      return State::prepare(frame);
    }
  }

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return Pattern::grid(
        Pattern::PlanePosition{projections::stereo_pattern_args(
            input.coords, params.pattern_freq)},
        source_kernel_params<GridSourceParams>(params), prepared);
  }
};

/**
 * @brief PrimitiveLattice source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame) accessors.
 */
template <typename State> struct PrimitiveLattice : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the lattice
   *        parameters.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ParamsProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        {
          State::params(frame).lattice_cell_scale
        } -> std::convertible_to<float>;
        {
          State::params(frame).lattice_shape_blend
        } -> std::convertible_to<float>;
        { State::params(frame).lattice_softness } -> std::convertible_to<float>;
        { State::params(frame).lattice_radius } -> std::convertible_to<float>;
      };

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame) {
    return Pattern::primitive_lattice(
        Pattern::PlanePosition{input.coords},
        source_kernel_params<LatticeSourceParams>(State::params(frame)));
  }
};

/**
 * @brief Tessellation source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), prepare(frame) accessors.
 */
template <typename State, Pattern::TessellationKind KindV>
struct Tessellation : PrepareFromState<State> {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the line
   *        parameters and rotation.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).cell_scale } -> std::convertible_to<float>;
        { State::params(frame).line_thickness } -> std::convertible_to<float>;
        { State::params(frame).line_softness } -> std::convertible_to<float>;
        { State::prepare(frame).angle_cos } -> std::convertible_to<float>;
        { State::prepare(frame).angle_sin } -> std::convertible_to<float>;
      };

  /// Phases `State::prepare` resolves for the frame.
  using Prepared = StatePrepared<State>;

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    return Pattern::tessellation(Pattern::PlanePosition{input.coords},
                                 State::params(frame), KindV, prepared);
  }
};

/**
 * @brief ProjectedNoise source policy.
 * @tparam State Provider with Binding and FrameState types and
 * noise(frame), noise_scale(frame), noise_time(frame),
 * noise_contrast(frame) accessors.
 */
template <typename State, math::NoiseBasis BasisV>
struct ProjectedNoise : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the noise
   *        generator, scale, time and contrast.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::noise_scale(frame) } -> std::same_as<float>;
        { State::noise_time(frame) } -> std::same_as<float>;
        { State::noise_contrast(frame) } -> std::same_as<float>;
      };

  /// Noise-space loop offset for the frame's noise time.
  using Prepared = math::Vector;

  /**
   * @brief Resolves the loop offset for the frame's noise time.
   * @param frame Frame state.
   * @return Offset added to the noise-space coordinate.
   */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return math::noise_projected_loop_offset(State::noise_time(frame));
  }

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    return Pattern::noise_contour(
        State::noise(frame), BasisV,
        math::noise_projected_coordinate(input.coords,
                                         State::noise_scale(frame), prepared),
        State::noise_contrast(frame));
  }
};

} // namespace Source

} // namespace Pullback
