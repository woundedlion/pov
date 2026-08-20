/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"
#include "render/pullback/fields.h"

/**
 * @file material.h
 * @brief Weight, transfer and coverage policies.
 */

namespace Pullback {

namespace Weight {

struct None : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(float field, const ProjectionProvenance &, const FrameState &) {
    return field;
  }
};

struct Projection : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(float field, const ProjectionProvenance &provenance,
        const FrameState &) {
    return field * provenance.value_weight;
  }
};

} // namespace Weight

namespace Transfer {

/** @brief Value placeholder for a chain whose value family has no fields. */
struct LinearValueParams {
  static constexpr std::array<Field<LinearValueParams>, 0> FIELDS{};
};
static_assert(field_ids_unique<LinearValueParams>());

/**
 * @brief Value parameters for the iso band
 *        (Pullback::Transfer::IsoContour).
 */
struct IsoValueParams {
  float iso_level = 0.5f;  /**< Source value the band is centered on. */
  float iso_width = 0.05f; /**< Half-width of the band's plateau. */

  static constexpr auto FIELDS = std::array{
      Field<IsoValueParams>{"iso-level", &IsoValueParams::iso_level,
                            "Iso Level", 0.0f, 1.0f, FieldCurve::LERP},
      Field<IsoValueParams>{"iso-width", &IsoValueParams::iso_width,
                            "Iso Width", 1.0f / 1024.0f, 1.0f,
                            FieldCurve::LOG_POSITIVE},
  };
};
static_assert(field_ids_unique<IsoValueParams>());

struct Linear : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static float apply(float value,
                                                    const FrameState &) {
    return value;
  }
};

struct Ridge : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static float apply(float value,
                                                    const FrameState &) {
    return unit_bell(value);
  }
};

/** @brief Shared iso-band kernel: a unit plateau of half-width @p width
    around @p level, falling to zero across a second half-width. */
__attribute__((always_inline)) inline float
iso_contour(float value, float level, float width) {
  const float distance = fabsf(value - level);
  return 1.0f - Detail::smooth_ramp(width, 2.0f * width, distance);
}

/** @brief Shared banding kernel: @p band_count cosine bands over the unit
    value, offset by @p band_phase. */
__attribute__((always_inline)) inline float
smooth_bands(float value, float band_count, float band_phase) {
  return 0.5f - 0.5f * cosf(TWO_PI_F * band_count * value + band_phase);
}

template <typename State> struct IsoContour : ApproximationDefaults {
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::iso_level(frame) } -> std::same_as<float>;
        { State::iso_width(frame) } -> std::same_as<float>;
      };

  template <typename FrameState>
  __attribute__((always_inline)) static float apply(float value,
                                                    const FrameState &frame) {
    return iso_contour(value, State::iso_level(frame), State::iso_width(frame));
  }
};

template <typename State> struct SmoothBands : ApproximationDefaults {
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::band_count(frame) } -> std::same_as<float>;
        { State::band_phase(frame) } -> std::same_as<float>;
      };

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
namespace ProjectedCoverage {

struct None : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(const ProjectionProvenance &, const FrameState &) {
    return 1.0f;
  }
};

struct Weight : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(const ProjectionProvenance &provenance, const FrameState &) {
    return provenance.value_weight;
  }
};

struct WeightSquared : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(const ProjectionProvenance &provenance, const FrameState &) {
    return provenance.value_weight * provenance.value_weight;
  }
};

/** @brief Shared edge-fade kernel; width 0 makes the edge a hard cut. */
__attribute__((always_inline)) inline float
edge_fade(const ProjectionProvenance &provenance, float width) {
  return width == 0.0f
             ? static_cast<float>(provenance.fade_edge_distance > 0.0f)
             : Detail::smooth_ramp(0.0f, width, provenance.fade_edge_distance);
}

template <typename State> struct EdgeFade : ApproximationDefaults {
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::edge_width(frame) } -> std::same_as<float>;
      };

  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(const ProjectionProvenance &provenance, const FrameState &frame) {
    return edge_fade(provenance, State::edge_width(frame));
  }
};

} // namespace ProjectedCoverage

namespace Coverage {

/** @brief Value parameters for the edge fade
    (Pullback::ProjectedCoverage::EdgeFade). */
struct EdgeValueParams {
  /** Fade band width in the projection's edge-distance units; 0 makes the edge
      a hard cut. */
  float edge_width = 0.1f;

  static constexpr auto FIELDS = std::array{
      Field<EdgeValueParams>{"edge-width", &EdgeValueParams::edge_width,
                             "Edge Width", 0.0f, 1.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<EdgeValueParams>());

/** @brief Shared cutout kernel: a smooth step through @p threshold with a
    half-width of @p width. */
__attribute__((always_inline)) inline float
value_cutout(float value, float threshold, float width) {
  return ::smooth_ramp(threshold - width, threshold + width, value);
}

/**
 * @brief Value-dependent coverage cut for Stage::Coverage.
 * @details Reads the current FIELD value and nothing else, so a chain may
 * legally place it before, between, or after transfers.
 */
template <typename State> struct ValueCutout : ApproximationDefaults {
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::cutout_threshold(frame) } -> std::same_as<float>;
        { State::cutout_width(frame) } -> std::same_as<float>;
      };

  template <typename FrameState>
  __attribute__((always_inline)) static float apply(float value,
                                                    const FrameState &frame) {
    return value_cutout(value, State::cutout_threshold(frame),
                        State::cutout_width(frame));
  }
};

} // namespace Coverage

} // namespace Pullback
