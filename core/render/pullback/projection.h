/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"
#include "render/pullback/fields.h"

/**
 * @file projection.h
 * @brief Sphere-to-plane projection policies.
 */

namespace Pullback {

namespace Projection {

/** @brief Projection and camera parameters, shared by every composed effect. */
struct ProjectionParams {
  float singularity_fade = 1.0f; /**< Radius of the projection's singularity
                                attenuation. */
  float spin_rate = 0.0f; /**< Per-frame spin of the projection frame about Y;
                                only read under `ANIMATED_PROJECTION`. */
  float wander = 0.0f;    /**< Fraction of the projection random-walk delta
                                absorbed each frame. */
  float camera_wander = 0.0f;    /**< Same, for the outer camera random walk. */
  float central_meridian = 0.0f; /**< Central meridian handed to projections
                                      that take one, in radians. */

  static constexpr auto FIELDS = std::array{
      Field<ProjectionParams>{"singularity-fade",
                              &ProjectionParams::singularity_fade,
                              "Singularity Fade", 1.0f, 20.0f, FieldCurve::LERP,
                              FieldGate::SINGULARITY_FADE},
      Field<ProjectionParams>{"projection-spin-speed",
                              &ProjectionParams::spin_rate,
                              "Projection Spin Speed", 0.0f, 0.05f,
                              FieldCurve::LERP, FieldGate::ANIMATED_PROJECTION},
      Field<ProjectionParams>{"projection-wander", &ProjectionParams::wander,
                              "Projection Wander", 0.0f, 1.0f, FieldCurve::LERP,
                              FieldGate::ANIMATED_PROJECTION},
      Field<ProjectionParams>{"camera-wander", &ProjectionParams::camera_wander,
                              "Camera Wander", 0.0f, 1.0f, FieldCurve::LERP},
      Field<ProjectionParams>{
          "central-meridian", &ProjectionParams::central_meridian,
          "Central Meridian", 0.0f, TWO_PI_F, FieldCurve::SHORTEST_PERIODIC,
          FieldGate::CENTRAL_MERIDIAN},
  };
};
static_assert(field_ids_unique<ProjectionParams>());

enum class GnomonicHemisphere : uint8_t { FOLDED, FRONT, BACK };

inline constexpr uint8_t FOLDED_FLAG = 1U << 0;
/** Render-space divisor floor that caps gnomonic coordinates near 1000. The
 * math primitive uses STEREO_EQUATOR_EPS only to avoid division by zero before
 * clamping to its much larger point-at-infinity sentinel. */
inline constexpr float GNOMONIC_AXIS_EPS = 1e-3f;

__attribute__((always_inline)) inline float
singularity_attenuation(float regular_distance_sq, float singular_distance_sq,
                        float singularity_fade) {
  const float pf = singularity_fade > 1e-3f ? singularity_fade : 1e-3f;
  const float scaled_distance_sq = pf * pf * regular_distance_sq;
  return scaled_distance_sq / (scaled_distance_sq + singular_distance_sq);
}

__attribute__((always_inline)) inline float
equirectangular_weight(const Vector &input, float singularity_fade) {
  return singularity_attenuation(input.x * input.x + input.z * input.z,
                                 input.y * input.y, singularity_fade);
}

__attribute__((always_inline)) inline float
equirectangular_weight(float latitude, float singularity_fade) {
  const float cos_latitude = cosf(latitude);
  const float sin_latitude = sinf(latitude);
  return singularity_attenuation(cos_latitude * cos_latitude,
                                 sin_latitude * sin_latitude, singularity_fade);
}

__attribute__((always_inline)) inline float
peirce_weight(const Vector &input, float meridian_cos, float meridian_sin,
              float singularity_fade) {
  const float rotated_x = input.x * meridian_cos + input.z * meridian_sin;
  const float rotated_z = input.z * meridian_cos - input.x * meridian_sin;
  const float singular_cosine =
      (fabsf(rotated_z) + fabsf(rotated_x)) * 0.7071067811865475f;
  float sin_distance_sq =
      std::max(0.0f, 1.0f - singular_cosine * singular_cosine);
  if (input.y < 0.0f) {
    const float fold_sine =
        fabsf(fabsf(rotated_z) - fabsf(rotated_x)) * 0.7071067811865475f;
    sin_distance_sq = std::min(sin_distance_sq, fold_sine * fold_sine);
  }
  return singularity_attenuation(sin_distance_sq,
                                 std::max(0.0f, 1.0f - sin_distance_sq),
                                 singularity_fade);
}

__attribute__((always_inline)) inline float
peirce_weight(const Vector &input, float central_meridian,
              float singularity_fade) {
  return peirce_weight(input, cosf(central_meridian), sinf(central_meridian),
                       singularity_fade);
}

__attribute__((always_inline)) inline ProjectionResult
stereographic(const Vector &input, float singularity_fade) {
  const Complex coords = stereo(input);
  return {coords,
          {0, 0, static_cast<uint8_t>(ProjectionBoundary::SINGULAR),
           std::max(0.0f, 1.0f - input.y),
           singularity_attenuation(std::max(0.0f, 1.0f - input.y),
                                   std::max(0.0f, 1.0f + input.y),
                                   singularity_fade),
           0}};
}

// Out of line under Emscripten, inlined on every other target.
#if defined(__EMSCRIPTEN__)
__attribute__((noinline))
#else
__attribute__((always_inline))
#endif
inline ProjectionResult folded_sinusoidal(const Vector &input,
                                          float central_meridian) {
  const Complex coords =
      projections::folded_sinusoidal(input, central_meridian);
  return {
      coords,
      {static_cast<uint8_t>(input.z < 0.0f), 0, 0, 1.0f, 1.0f, FOLDED_FLAG}};
}

__attribute__((always_inline)) inline ProjectionResult
equirectangular(const Vector &input, float central_meridian,
                float singularity_fade) {
  const Complex coords = projections::equirectangular(input, central_meridian);
  return {coords,
          {0, 0, static_cast<uint8_t>(ProjectionBoundary::CUT),
           PI_F - fabsf(coords.re),
           equirectangular_weight(input, singularity_fade), 0}};
}

__attribute__((always_inline)) inline ProjectionResult
gnomonic(const Vector &input, float singularity_fade,
         GnomonicHemisphere hemisphere) {
  float y = input.y;
  if (fabsf(y) < GNOMONIC_AXIS_EPS)
    y = y < 0.0f ? -GNOMONIC_AXIS_EPS : GNOMONIC_AXIS_EPS;
  const Complex coords(input.x / y, input.z / y);
  const bool in_domain =
      hemisphere == GnomonicHemisphere::FOLDED ||
      (hemisphere == GnomonicHemisphere::FRONT ? input.y >= 0.0f
                                               : input.y < 0.0f);
  return {
      coords,
      {static_cast<uint8_t>(input.y < 0.0f),
       static_cast<uint8_t>(input.y < 0.0f),
       static_cast<uint8_t>(static_cast<uint8_t>(ProjectionBoundary::CUT) |
                            static_cast<uint8_t>(ProjectionBoundary::SINGULAR)),
       fabsf(input.y),
       singularity_attenuation(input.y * input.y,
                               input.x * input.x + input.z * input.z,
                               singularity_fade),
       0, 0, 0, in_domain ? 1.0f : 0.0f}};
}

__attribute__((always_inline)) inline ProjectionResult
from_kernel(const projections::ProjectionKernelResult &result,
            float coordinate_scale, float value_weight = 1.0f) {
  return {{result.coords.re * coordinate_scale,
           result.coords.im * coordinate_scale},
          {result.region_id, result.component_id, result.boundary_flags,
           result.fade_edge_distance * fabsf(coordinate_scale), value_weight,
           result.flags, result.traits, result.edge_class}};
}

__attribute__((always_inline)) inline ProjectionResult
bonne(const Vector &input, float central_meridian, float standard_parallel,
      float coordinate_scale) {
  return from_kernel(
      projections::bonne_projection(input, central_meridian, standard_parallel),
      coordinate_scale);
}

__attribute__((always_inline)) inline ProjectionResult
peirce(const Vector &input, float central_meridian, uint8_t layout,
       float layout_scroll, bool edge_distance_required, float coordinate_scale,
       float singularity_fade, float meridian_cos, float meridian_sin) {
  return from_kernel(
      projections::peirce_projection(input, central_meridian, layout,
                                     layout_scroll, edge_distance_required),
      coordinate_scale,
      peirce_weight(input, meridian_cos, meridian_sin, singularity_fade));
}

__attribute__((always_inline)) inline ProjectionResult
peirce(const Vector &input, float central_meridian, uint8_t layout,
       float layout_scroll, bool edge_distance_required, float coordinate_scale,
       float singularity_fade) {
  return peirce(input, central_meridian, layout, layout_scroll,
                edge_distance_required, coordinate_scale, singularity_fade,
                cosf(central_meridian), sinf(central_meridian));
}

__attribute__((always_inline)) inline ProjectionResult
peirce_fast_square(const Vector &input, float coordinate_scale,
                   float singularity_fade) {
  return from_kernel(projections::peirce_projection_fast_square(input),
                     coordinate_scale,
                     peirce_weight(input, 0.0f, singularity_fade));
}

template <typename State, typename Binding>
concept FrameProvider = Detail::ProviderFor<State, Binding> &&
                        requires(const typename Binding::FrameState &frame) {
                          {
                            State::conjugate(frame)
                          } -> std::same_as<const Quaternion &>;
                        };

__attribute__((always_inline)) inline ProjectionResult
airocean(const Vector &input, float central_meridian, bool horizontal,
         bool edge_distance_required, float coordinate_scale) {
  return from_kernel(projections::airocean_projection(input, central_meridian,
                                                      horizontal,
                                                      edge_distance_required),
                     coordinate_scale);
}

/** @brief Bonne pseudoconical equal-area projection; `North` picks the sign of
    the standard parallel, and so the hemisphere the cone opens toward. */
template <typename State, bool North> struct Bonne : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::standard_parallel(frame) } -> std::same_as<float>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    const float hemisphere = North ? 1.0f : -1.0f;
    return bonne(input, State::central_meridian(frame),
                 hemisphere * State::standard_parallel(frame),
                 State::coordinate_scale(frame));
  }
};

/** @brief Stereographic projection: conformal, with one singular pole the
    singularity fade attenuates. */
template <typename State> struct Stereographic : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::singularity_fade(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return stereographic(input, State::singularity_fade(frame));
  }
};

/** @brief Sinusoidal projection with the azimuth folded about the central
    meridian: both hemispheres share one image, and there is no singular locus
    to attenuate. */
template <typename State> struct FoldedSinusoidal : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return folded_sinusoidal(input, State::central_meridian(frame));
  }
};

/** @brief Equirectangular projection: cut at the antimeridian, with both
    poles attenuated by the singularity fade. */
template <typename State> struct Equirectangular : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::singularity_fade(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return equirectangular(input, State::central_meridian(frame),
                           State::singularity_fade(frame));
  }
};

/** @brief Gnomonic projection about the Y axis, singular on the y = 0 great
    circle; `Hemisphere` folds the two halves together or keeps one. */
template <typename State, GnomonicHemisphere Hemisphere>
struct Gnomonic : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::singularity_fade(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return gnomonic(input, State::singularity_fade(frame), Hemisphere);
  }
};

/** @brief Peirce quincuncial projection, conformal but for four singularities;
    `Layout` picks diamond, square or strip tiling and `EdgeDistanceRequired`
    makes the kernel compute edge distance unconditionally. */
template <typename State, uint8_t Layout, bool EdgeDistanceRequired>
struct Peirce : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  static constexpr bool EDGE_DISTANCE_UNCONDITIONAL = EdgeDistanceRequired;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::layout_scroll(frame) } -> std::same_as<float>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
        { State::singularity_fade(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return peirce(input, State::central_meridian(frame), Layout,
                  State::layout_scroll(frame), EdgeDistanceRequired,
                  State::coordinate_scale(frame),
                  State::singularity_fade(frame));
  }
};

/** @brief Approximation bounds of the fast square Peirce path, shared by the
    template policy and the chain operator's descriptor. */
inline constexpr std::array<ApproximationMetric, 3> PEIRCE_FAST_SQUARE_METRICS{{
    {ApproximationDomain::PROJECTED_COORDINATE,
     ApproximationAggregation::MAXIMUM, 1.2e-3f, "plane units"},
    {ApproximationDomain::PROJECTED_EDGE_DISTANCE,
     ApproximationAggregation::MAXIMUM, 2e-4f, "plane units"},
    {ApproximationDomain::FRAMEBUFFER, ApproximationAggregation::MAXIMUM,
     256.0f, "channel code"},
}};

/** @brief Approximate square-layout Peirce projection; the provider must pin
    the central meridian to zero. */
template <typename State> struct PeirceFastSquare : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  static constexpr bool EDGE_DISTANCE_UNCONDITIONAL = true;
  static constexpr bool APPROXIMATE = true;
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::PEIRCE_FAST_SQUARE;
  static constexpr std::array<ApproximationMetric, 3> METRICS =
      PEIRCE_FAST_SQUARE_METRICS;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::ZERO_CENTRAL_MERIDIAN } -> std::convertible_to<bool>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
        { State::singularity_fade(frame) } -> std::same_as<float>;
      } && State::ZERO_CENTRAL_MERIDIAN;

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return peirce_fast_square(input, State::coordinate_scale(frame),
                              State::singularity_fade(frame));
  }

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }
};

/**
 * @brief Square-layout Peirce projection taking the approximate path only at a
 *        zero central meridian.
 * @details The approximation oracle and metrics are inherited from
 * PeirceFastSquare and so cover the whole policy; off a zero central meridian
 * it runs the exact quincuncial kernel and those bounds are slack.
 */
template <typename State> struct PeirceSquare : PeirceFastSquare<State> {
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::layout_scroll(frame) } -> std::same_as<float>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
        { State::singularity_fade(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    if (State::central_meridian(frame) == 0.0f)
      return peirce_fast_square(input, State::coordinate_scale(frame),
                                State::singularity_fade(frame));
    return peirce(
        input, State::central_meridian(frame), 1, State::layout_scroll(frame),
        true, State::coordinate_scale(frame), State::singularity_fade(frame));
  }
};

/** @brief Airocean icosahedral net; `Horizontal` turns the finished net a
    quarter turn and `EdgeDistanceRequired` makes the kernel compute the
    per-edge cut distances unconditionally. */
template <typename State, bool Horizontal, bool EdgeDistanceRequired>
struct Airocean : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  static constexpr bool EDGE_DISTANCE_UNCONDITIONAL = EdgeDistanceRequired;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return airocean(input, State::central_meridian(frame), Horizontal,
                    EdgeDistanceRequired, State::coordinate_scale(frame));
  }
};

} // namespace Projection

} // namespace Pullback
