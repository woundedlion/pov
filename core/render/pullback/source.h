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

/**
 * @brief Source parameters for the coupled sine grid
 *        (Pullback::Source::Grid).
 */
struct GridSourceParams : Pattern::GridSourceParams {
  /// translation_lattice() is axis-aligned for every parameter value.
  static constexpr bool RECTANGULAR_TRANSLATION_LATTICE = true;
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
  };
};
static_assert(field_ids_unique<GridSourceParams>());
static_assert(field_defaults_in_range<GridSourceParams>());

/** @brief Square translation symmetry of the coupled grid.
    @details Holds wherever the frequency-scaled coordinates stay inside
    `projections::STEREO_PATTERN_ARG_LIMIT`, beyond which the source clamps.
    @param params Grid controls.
    @return Square lattice of period `2pi / pattern_freq`. */
constexpr TranslationLattice
translation_lattice(const GridSourceParams &params) {
  return square_lattice(math::TWO_PI_F / params.pattern_freq);
}

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

/** @brief Source parameters for parallel plane bands. */
struct BandsSourceParams : Pattern::BandsSourceParams {
  /** @brief Numeric source controls. */
  static constexpr auto FIELDS = std::array{
      Field<BandsSourceParams>{"pattern-freq", &BandsSourceParams::pattern_freq,
                               "Pattern Freq", 0.01f, 64.0f, FieldCurve::LERP},
      Field<BandsSourceParams>{"angle", &BandsSourceParams::angle, "Angle",
                               0.0f, math::TWO_PI_F,
                               FieldCurve::SHORTEST_PERIODIC},
  };
};
static_assert(field_ids_unique<BandsSourceParams>());
static_assert(field_defaults_in_range<BandsSourceParams>());

/** @brief Translation symmetry of parallel bands within the unclamped source domain.
 * @param params Band frequency and orientation.
 * @return A rotated square lattice with one band-normal period. */
inline TranslationLattice translation_lattice(const BandsSourceParams &params) {
  const float PERIOD = math::TWO_PI_F / params.pattern_freq;
  const float COSINE = cosf(params.angle), SINE = sinf(params.angle);
  return {{PERIOD * COSINE, PERIOD * SINE}, {-PERIOD * SINE, PERIOD * COSINE}};
}

/**
 * @brief Source parameters for the stationary spiral field
 *        (Pullback::Source::Spiral).
 */
struct SpiralSourceParams : Pattern::SpiralSourceParams {
  /** @brief Per-parameter id, member, label, range and curve. */
  static constexpr auto FIELDS = std::array{
      Field<SpiralSourceParams>{"pattern-freq",
                                &SpiralSourceParams::pattern_freq,
                                "Pattern Freq", 0.1f, 20.0f, FieldCurve::LERP},
      Field<SpiralSourceParams>{"arms", &SpiralSourceParams::arms,
                                "Spiral Arms", 0.0f, 12.0f, FieldCurve::SNAP},
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

/** @brief Fill topology gate of the cell-primitive controls. */
inline constexpr TopologyGate TESSELLATION_CELLS_GATE{
    "fill", live_values(Pattern::TessellationFill::CELLS)};
/** @brief Fill topology gate of the edge-band controls. */
inline constexpr TopologyGate TESSELLATION_EDGES_GATE{
    "fill", live_values(Pattern::TessellationFill::EDGES)};

/** @brief Field table of the tessellation controls @p Owner carries.
 * @details The cell and edge controls are gated on the fill topology only when
 * @p Owner carries both.
 * @tparam Owner Tessellation parameter family.
 * @tparam Members Complete base of @p Owner declaring its controls.
 * @return Cell scale, then the cell-primitive and edge-band fields. */
template <typename Owner, typename Members>
consteval auto tessellation_fields() {
  constexpr bool CELLS = requires(Members members) { members.cell_radius; };
  constexpr bool EDGES = requires(Members members) { members.line_thickness; };
  static_assert(CELLS || EDGES);
  const Field<Owner> scale{"cell-scale", &Owner::cell_scale,
                           "Cell Scale", 1.0f / 64.0f,
                           8.0f,         FieldCurve::LOG_POSITIVE};
  const auto cells = [](TopologyGate gate) {
    if constexpr (CELLS)
      return std::array{
          Field<Owner>{"cell-shape", &Owner::cell_shape, "Cell Shape", 0.0f,
                       1.0f, FieldCurve::LERP, FieldGate::ALWAYS, gate},
          Field<Owner>{"cell-softness", &Owner::cell_softness, "Cell Softness",
                       1.0f / 1024.0f, 1.0f, FieldCurve::LOG_POSITIVE,
                       FieldGate::ALWAYS, gate},
          Field<Owner>{"cell-radius", &Owner::cell_radius, "Cell Radius",
                       1.0f / 64.0f, 0.49f, FieldCurve::LERP, FieldGate::ALWAYS,
                       gate},
      };
    else
      return std::array<Field<Owner>, 0>{};
  };
  const auto edges = [](TopologyGate gate) {
    if constexpr (EDGES)
      return std::array{
          Field<Owner>{"line-thickness", &Owner::line_thickness,
                       "Line Thickness", 1.0f / 1024.0f, 0.25f,
                       FieldCurve::LOG_POSITIVE, FieldGate::ALWAYS, gate},
          Field<Owner>{"line-softness", &Owner::line_softness, "Line Softness",
                       1.0f / 1024.0f, 0.25f, FieldCurve::LOG_POSITIVE,
                       FieldGate::ALWAYS, gate},
      };
    else
      return std::array<Field<Owner>, 0>{};
  };
  constexpr bool GATED = CELLS && EDGES;
  const auto primitive =
      cells(GATED ? TESSELLATION_CELLS_GATE : TopologyGate{});
  const auto band = edges(GATED ? TESSELLATION_EDGES_GATE : TopologyGate{});
  std::array<Field<Owner>, 1 + primitive.size() + band.size()> out{};
  out[0] = scale;
  for (size_t index = 0; index < primitive.size(); ++index)
    out[1 + index] = primitive[index];
  for (size_t index = 0; index < band.size(); ++index)
    out[1 + primitive.size() + index] = band[index];
  return out;
}

/** @brief Translation symmetry of a tessellation.
    @param kind Cell shape.
    @param cell_scale Lattice periods per plane unit.
    @return Square lattice for square cells, else the hexagonal lattice, both
    of period `1 / cell_scale`. */
constexpr TranslationLattice
tessellation_lattice(Pattern::TessellationKind kind, float cell_scale) {
  const float period = 1.0f / cell_scale;
  return kind == Pattern::TessellationKind::SQUARE ? square_lattice(period)
                                                   : hexagonal_lattice(period);
}

/** @brief Source parameters for the static tessellation with both fills'
    controls, gated on a `fill` topology the deriving family declares. */
struct TessellationSourceParams : Pattern::TessellationSourceParams {
  /** @brief Per-parameter id, member, label, range and curve. */
  static constexpr auto FIELDS =
      tessellation_fields<TessellationSourceParams,
                          Pattern::TessellationSourceParams>();
};
static_assert(field_ids_unique<TessellationSourceParams>());
static_assert(field_defaults_in_range<TessellationSourceParams>());

/** @brief Pattern controls one tessellation fill reads.
 * @tparam Fill Cell primitives or edge band. */
template <Pattern::TessellationFill Fill>
using TessellationControls =
    std::conditional_t<Fill == Pattern::TessellationFill::CELLS,
                       Pattern::TessellationCellParams,
                       Pattern::TessellationEdgeParams>;

/**
 * @brief Tessellation source controls with a compile-time kind and fill
 *        (Pullback::Source::Tessellation).
 * @tparam Kind Cell shape.
 * @tparam Fill Cell primitives or edge band; selects the carried controls.
 */
template <Pattern::TessellationKind Kind, Pattern::TessellationFill Fill>
struct FixedTessellationParams : TessellationControls<Fill> {
  static constexpr Pattern::TessellationKind KIND = Kind; ///< Cell shape.
  static constexpr Pattern::TessellationFill FILL = Fill; ///< Fill.
  /// translation_lattice() is axis-aligned for every parameter value.
  static constexpr bool RECTANGULAR_TRANSLATION_LATTICE =
      Kind == Pattern::TessellationKind::SQUARE;
  /** @brief Per-parameter id, member, label, range and curve. */
  static constexpr auto FIELDS =
      tessellation_fields<FixedTessellationParams,
                          TessellationControls<Fill>>();
};

/** @brief Square-cell tessellation controls. */
using SquareCellsParams =
    FixedTessellationParams<Pattern::TessellationKind::SQUARE,
                            Pattern::TessellationFill::CELLS>;
static_assert(field_ids_unique<SquareCellsParams>());
static_assert(field_defaults_in_range<SquareCellsParams>());

/** @brief Translation symmetry of a fixed tessellation.
    @tparam Kind Cell shape.
    @tparam Fill Fill.
    @param params Tessellation controls.
    @return `tessellation_lattice` of @p Kind. */
template <Pattern::TessellationKind Kind, Pattern::TessellationFill Fill>
constexpr TranslationLattice
translation_lattice(const FixedTessellationParams<Kind, Fill> &params) {
  return tessellation_lattice(Kind, params.cell_scale);
}

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

/** @brief Parallel bands in the projection plane.
 * @tparam State Provider with Binding, FrameState and params(frame). */
template <typename State> struct Bands : ApproximationDefaults {
  using Binding = typename State::Binding;       ///< Provider binding.
  using FrameState = typename State::FrameState; ///< Provider frame state.
  using Prepared =
      Pattern::PlaneBands::Prepared; ///< Frequency and orientation.

  /** @brief Whether the provider supplies band controls.
   * @tparam CandidateBinding Binding being checked. */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ParamsProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::params(frame).angle } -> std::convertible_to<float>;
      };

  /** @brief Prepares band orientation.
   * @param frame Frame state.
   * @return Immutable band inputs. */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    const auto &params = State::params(frame);
    return Pattern::PlaneBands::prepare({params.pattern_freq, params.angle},
                                        {});
  }
  /** @brief Samples the projected plane.
   * @param input Plane carrier.
   * @param frame Frame state.
   * @param prepared Immutable band inputs.
   * @return Signed field value in [-1, 1]. */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    (void)frame;
    return Pattern::PlaneBands::sample({input.coords}, prepared);
  }
};

/**
 * @brief Stationary spiral source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame) accessors.
 */
template <typename State> struct Spiral : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the pattern
   *        frequency and arm count.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ParamsProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::params(frame).arms } -> std::convertible_to<float>;
      };

  /** @brief The whole arm count. */
  struct Prepared {
    float arms; ///< `arms` rounded to whole.
  };

  /** @brief Rounds the arm count once per frame.
   * @param frame Frame state.
   * @return Whole arm count. */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return {std::round(State::params(frame).arms)};
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
    return Pattern::spiral(
        Pattern::PlanePosition{projections::stereo_pattern_args(
            input.coords, State::params(frame).pattern_freq)},
        prepared.arms);
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
   *        parameters and coupling phase.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::params(frame).pattern_mix } -> std::convertible_to<float>;
        { State::params(frame).complexity } -> std::convertible_to<float>;
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
    return Pattern::grid(
        Pattern::PlanePosition{projections::stereo_pattern_args(
            input.coords, params.pattern_freq)},
        source_kernel_params<GridSourceParams>(params), prepared);
  }
};

/**
 * @brief Tessellation source policy.
 * @tparam State Provider with Binding and FrameState types and a
 * params(frame) accessor.
 * @tparam KindV Cell shape.
 * @tparam FillV Cell primitives or edge band.
 */
template <typename State, Pattern::TessellationKind KindV,
          Pattern::TessellationFill FillV>
struct Tessellation : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the controls
   *        `FillV` reads.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ParamsProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).cell_scale } -> std::convertible_to<float>;
      } &&
      (FillV == Pattern::TessellationFill::CELLS ? requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).cell_shape } -> std::convertible_to<float>;
        { State::params(frame).cell_softness } -> std::convertible_to<float>;
        { State::params(frame).cell_radius } -> std::convertible_to<float>;
      } : requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).line_thickness } -> std::convertible_to<float>;
        { State::params(frame).line_softness } -> std::convertible_to<float>;
      });

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame) {
    return Pattern::tessellation<KindV, FillV>(
        Pattern::PlanePosition{input.coords},
        source_kernel_params<TessellationControls<FillV>>(
            State::params(frame)));
  }
};

/**
 * @brief ProjectedNoise source policy.
 * @tparam State Provider with Binding and FrameState types and
 * noise(frame), noise_scale(frame), noise_time(frame),
 * noise_contrast(frame) accessors.
 * @tparam BasisV Noise basis.
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
