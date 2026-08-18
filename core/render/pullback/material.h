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

struct None : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(float field, const ProjectionSample &, const FrameState &) {
    return field;
  }
};

struct Projection : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(float field, const ProjectionSample &projected, const FrameState &) {
    return field * projected.value_weight;
  }
};

} // namespace Weight

namespace Transfer {

/** @brief Value placeholder for a material stage that takes no parameters. */
struct LinearValueParams {
  static constexpr std::array<Field<LinearValueParams>, 0> FIELDS{};
};

/**
 * @brief Value parameters for the iso band
 *        (Pullback::Transfer::IsoContour).
 */
struct IsoValueParams {
  float iso_level = 0.5f;  /**< Source value the band is centered on. */
  float iso_width = 0.05f; /**< Half-width of the band's plateau. */

  static constexpr auto FIELDS = std::array{
      Field<IsoValueParams>{&IsoValueParams::iso_level, "Iso Level", 0.0f, 1.0f,
                            FieldCurve::LERP},
      Field<IsoValueParams>{&IsoValueParams::iso_width, "Iso Width",
                            1.0f / 1024.0f, 1.0f, FieldCurve::LOG_POSITIVE},
  };
};

struct Linear : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static float apply(float value,
                                                    const FrameState &) {
    return value;
  }
};

struct Ridge : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static float apply(float value,
                                                    const FrameState &) {
    return unit_bell(value);
  }
};

template <typename State> struct IsoContour : ExactPolicy {
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
    const float width = State::iso_width(frame);
    const float distance = fabsf(value - State::iso_level(frame));
    return 1.0f - Detail::smooth_ramp(width, 2.0f * width, distance);
  }
};

template <typename State> struct SmoothBands : ExactPolicy {
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
    return 0.5f - 0.5f * cosf(TWO_PI_F * State::band_count(frame) * value +
                              State::band_phase(frame));
  }
};

} // namespace Transfer

namespace Coverage {

/** @brief Value parameters for the edge fade (Pullback::Coverage::EdgeFade). */
struct EdgeValueParams {
  /** Fade band width in the projection's edge-distance units; 0 makes the edge
      a hard cut. */
  float edge_width = 0.1f;

  static constexpr auto FIELDS = std::array{
      Field<EdgeValueParams>{&EdgeValueParams::edge_width, "Edge Width", 0.0f,
                             1.0f, FieldCurve::LERP},
  };
};

struct Opaque : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(float, const ProjectionSample &, const FrameState &) {
    return 1.0f;
  }
};

struct Projection : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(float, const ProjectionSample &projected, const FrameState &) {
    return projected.value_weight;
  }
};

struct ProjectionSquared : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(float, const ProjectionSample &projected, const FrameState &) {
    return projected.value_weight * projected.value_weight;
  }
};

template <typename State> struct ValueCutout : ExactPolicy {
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::cutout_threshold(frame) } -> std::same_as<float>;
        { State::cutout_width(frame) } -> std::same_as<float>;
      };

  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(float value, const ProjectionSample &, const FrameState &frame) {
    const float threshold = State::cutout_threshold(frame);
    const float width = State::cutout_width(frame);
    return ::smooth_ramp(threshold - width, threshold + width, value);
  }
};

template <typename State> struct EdgeFade : ExactPolicy {
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::edge_width(frame) } -> std::same_as<float>;
      };

  template <typename FrameState>
  __attribute__((always_inline)) static float
  apply(float, const ProjectionSample &projected, const FrameState &frame) {
    const float width = State::edge_width(frame);
    return width == 0.0f
               ? static_cast<float>(projected.fade_edge_distance > 0.0f)
               : Detail::smooth_ramp(0.0f, width, projected.fade_edge_distance);
  }
};

} // namespace Coverage

} // namespace Pullback
