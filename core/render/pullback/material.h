/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"
#include "render/pullback/fields.h"
#include "math/3dmath.h"

/**
 * @file material.h
 * @brief Weight, transfer and coverage policies.
 */

namespace Pullback {

/** @brief Projection-coverage modes consumed by the Sample crossing. */
enum class ProjectionCoverageMode : uint8_t {
  NONE = 0,
  WEIGHT = 1,
  WEIGHT_SQUARED = 2,
  EDGE_FADE = 3
};

namespace Weight {

/** @brief Weight policy that leaves the field unweighted. */
struct None : ApproximationDefaults {
  /**
   * @brief Passes the field through.
   * @tparam FrameState Frame state type; unused.
   * @param field Field value.
   * @return `field`.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(float field, const ProjectionProvenance &, const FrameState &) {
    return field;
  }
};

/** @brief Weight policy scaling the field by the projection weight. */
struct Projection : ApproximationDefaults {
  /**
   * @brief Scales the field by the projection's value weight.
   * @tparam FrameState Frame state type; unused.
   * @param field Field value.
   * @param provenance Projection provenance of the sample.
   * @return `field * provenance.value_weight`.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(float field, const ProjectionProvenance &provenance,
        const FrameState &) {
    return field * provenance.value_weight;
  }
};

} // namespace Weight

namespace Transfer {

/** @brief Value placeholder for a chain with no value-stage parameters. */
struct NoValueParams {
  /// Empty parameter registry.
  static constexpr std::array<Field<NoValueParams>, 0> FIELDS{};
};
static_assert(field_ids_unique<NoValueParams>());
static_assert(field_defaults_in_range<NoValueParams>());

/**
 * @brief Value parameters for the iso band
 *        (Pullback::Transfer::IsoContour).
 */
struct IsoValueParams {
  float iso_level = 0.5f;  /**< Source value the band is centered on. */
  float iso_width = 0.05f; /**< Half-width of the band's plateau. */

  /// Parameter registry: id, member, label, range and curve per field.
  static constexpr auto FIELDS = std::array{
      Field<IsoValueParams>{"iso-level", &IsoValueParams::iso_level,
                            "Iso Level", 0.0f, 1.0f, FieldCurve::LERP},
      Field<IsoValueParams>{"iso-width", &IsoValueParams::iso_width,
                            "Iso Width", 1.0f / 1024.0f, 1.0f,
                            FieldCurve::LOG_POSITIVE},
  };
};
static_assert(field_ids_unique<IsoValueParams>());
static_assert(field_defaults_in_range<IsoValueParams>());

/** @brief Transfer peaking at mid-value: `math::unit_bell`. */
struct Ridge : ApproximationDefaults, TransferRole {
  /**
   * @brief Applies the ridge transfer.
   * @tparam FrameState Frame state type; unused.
   * @param value Field value in [0, 1].
   * @return `math::unit_bell(value)`.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static float apply(float value,
                                                    const FrameState &) {
    return math::unit_bell(value);
  }
};

/** @brief Shared iso-band kernel: a unit plateau of half-width @p width
    around @p level, falling to zero across a second half-width. */
__attribute__((always_inline)) inline float
iso_contour(float value, float level, float width) {
  const float distance = fabsf(value - level);
  return 1.0f - Detail::smooth_ramp_or_step(width, 2.0f * width, distance);
}

/** @brief Shared banding kernel: @p band_count cosine bands over the unit
    value, offset by @p band_phase. */
__attribute__((always_inline)) inline float
smooth_bands(float value, float band_count, float band_phase) {
  return 0.5f - 0.5f * math::fast_cosf(math::TWO_PI_F * band_count * value +
                                       band_phase);
}

/**
 * @brief IsoContour material policy.
 * @tparam State Provider with Binding and FrameState types and
 * iso_level(frame), iso_width(frame) accessors.
 */
template <typename State>
struct IsoContour : ApproximationDefaults, TransferRole {
  /**
   * @brief Whether `State` is a provider for `Binding` with every accessor
   * this policy reads.
   * @tparam Binding Chain binding to check against.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::iso_level(frame) } -> std::same_as<float>;
        { State::iso_width(frame) } -> std::same_as<float>;
      };

  /**
   * @brief Applies this frame's iso band.
   * @tparam FrameState Frame state of the provider's binding.
   * @param value Field value.
   * @param frame Frame state the level and width are read from.
   * @return iso_contour() of `value`.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static float apply(float value,
                                                    const FrameState &frame) {
    return iso_contour(value, State::iso_level(frame), State::iso_width(frame));
  }
};

/**
 * @brief SmoothBands material policy.
 * @tparam State Provider with Binding and FrameState types and
 * band_count(frame), band_phase(frame) accessors.
 */
template <typename State>
struct SmoothBands : ApproximationDefaults, TransferRole {
  /**
   * @brief Whether `State` is a provider for `Binding` with every accessor
   * this policy reads.
   * @tparam Binding Chain binding to check against.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::band_count(frame) } -> std::same_as<float>;
        { State::band_phase(frame) } -> std::same_as<float>;
      };

  /**
   * @brief Applies this frame's cosine bands.
   * @tparam FrameState Frame state of the provider's binding.
   * @param value Field value.
   * @param frame Frame state the band count and phase are read from.
   * @return smooth_bands() of `value`, in [0, 1].
   */
  template <typename FrameState>
  __attribute__((always_inline)) static float apply(float value,
                                                    const FrameState &frame) {
    return smooth_bands(value, State::band_count(frame),
                        State::band_phase(frame));
  }
};

} // namespace Transfer

/**
 * @brief Coverage vocabulary of the Sample crossing: mutually exclusive modes
 *        over the projection provenance, consumed exactly once per chain.
 */
namespace ProjectionCoverage {

/** @brief Full coverage regardless of provenance. */
struct None : ApproximationDefaults {
  /**
   * @brief Full coverage.
   * @tparam FrameState Frame state type; unused.
   * @return 1.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(const ProjectionProvenance &, const FrameState &) {
    return 1.0f;
  }
};

/** @brief Coverage equal to the projection weight. */
struct Weight : ApproximationDefaults {
  /**
   * @brief Coverage from the projection weight.
   * @tparam FrameState Frame state type; unused.
   * @param provenance Projection provenance of the sample.
   * @return `provenance.value_weight`.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(const ProjectionProvenance &provenance, const FrameState &) {
    return provenance.value_weight;
  }
};

/** @brief Coverage equal to the squared projection weight. */
struct WeightSquared : ApproximationDefaults {
  /**
   * @brief Coverage from the squared projection weight.
   * @tparam FrameState Frame state type; unused.
   * @param provenance Projection provenance of the sample.
   * @return `provenance.value_weight` squared.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(const ProjectionProvenance &provenance, const FrameState &) {
    return provenance.value_weight * provenance.value_weight;
  }
};

using Detail::edge_fade;

/**
 * @brief EdgeFade material policy.
 * @tparam State Provider with Binding and FrameState types and
 * edge_width(frame) accessors.
 */
template <typename State> struct EdgeFade : ApproximationDefaults {
  /**
   * @brief Whether `State` is a provider for `Binding` with every accessor
   * this policy reads.
   * @tparam Binding Chain binding to check against.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::edge_width(frame) } -> std::same_as<float>;
      };

  /**
   * @brief Coverage fading toward the projection's fade-eligible edge.
   * @tparam FrameState Frame state of the provider's binding.
   * @param provenance Projection provenance of the sample.
   * @param frame Frame state the edge width is read from.
   * @return `edge_fade(provenance, edge_width)`.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(const ProjectionProvenance &provenance, const FrameState &frame) {
    return edge_fade(provenance, State::edge_width(frame));
  }
};

/** @brief Parameters for projection edge-fade coverage. */
struct EdgeValueParams {
  /** Fade band width in the projection's edge-distance units; 0 makes the edge
      a hard cut. */
  float edge_width = 0.1f;

  /// Parameter registry: id, member, label, range and curve per field.
  static constexpr auto FIELDS = std::array{
      edge_width_field(&EdgeValueParams::edge_width),
  };
};
static_assert(field_ids_unique<EdgeValueParams>());
static_assert(field_defaults_in_range<EdgeValueParams>());

} // namespace ProjectionCoverage

namespace ValueCoverage {

/** @brief Parameters for a value-dependent coverage cut. */
struct CutoutValueParams {
  float cutout_threshold = 0.5f; /**< Value the cut steps through. */
  float cutout_softness = 0.05f; /**< Half-width of the step's ramp. */

  /// Parameter registry: id, member, label, range and curve per field.
  static constexpr auto FIELDS = std::array{
      Field<CutoutValueParams>{
          "cutout-threshold", &CutoutValueParams::cutout_threshold,
          "Cutout Threshold", 0.0f, 1.0f, FieldCurve::LERP},
      Field<CutoutValueParams>{
          "cutout-softness", &CutoutValueParams::cutout_softness,
          "Cutout Softness", 1.0f / 1024.0f, 0.5f, FieldCurve::LOG_POSITIVE},
  };
};
static_assert(field_ids_unique<CutoutValueParams>());
static_assert(field_defaults_in_range<CutoutValueParams>());

/** @brief Shared cutout kernel: a smooth step through @p threshold with a
    half-width of @p width. */
__attribute__((always_inline)) inline float
value_cutout(float value, float threshold, float width) {
  return Detail::smooth_ramp_or_step(threshold - width, threshold + width,
                                     value);
}

/**
 * @brief Value-dependent coverage cut for Stage::ApplyCoverage.
 * @details Reads the current FIELD value and nothing else, so a chain may
 * legally place it before, between, or after transfers.
 * @tparam State Provider with Binding and FrameState types and
 * cutout_threshold(frame), cutout_softness(frame) accessors.
 */
template <typename State>
struct ValueCutout : ApproximationDefaults, CoverageRole {
  /**
   * @brief Whether `State` is a provider for `Binding` with every accessor
   * this policy reads.
   * @tparam Binding Chain binding to check against.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::cutout_threshold(frame) } -> std::same_as<float>;
        { State::cutout_softness(frame) } -> std::same_as<float>;
      };

  /**
   * @brief Coverage factor of this frame's cutout.
   * @tparam FrameState Frame state of the provider's binding.
   * @param value Field value.
   * @param frame Frame state the threshold and softness are read from.
   * @return value_cutout() of `value`, in [0, 1].
   */
  template <typename FrameState>
  __attribute__((always_inline)) static float apply(float value,
                                                    const FrameState &frame) {
    return value_cutout(value, State::cutout_threshold(frame),
                        State::cutout_softness(frame));
  }
};

} // namespace ValueCoverage

} // namespace Pullback
