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

/** @brief Canonical base orientation of a projection.
    @return Rotation taking -Z onto -Y. */
inline math::Quaternion projection_base_orientation() {
  return math::make_rotation(math::Vector(0, 0, -1), math::Vector(0, -1, 0));
}

namespace Projection {

/** @brief Projection and camera parameters. */
struct ProjectionParams {
  float singularity_fade = 1.0f; /**< Sharpness of the singularity attenuation:
                                      1 fades to the regular locus; 20 confines
                                      the fade to a narrow cap. */
  float spin_rate = 0.0f; /**< Per-frame spin of the projection frame about Y;
                                only read under `ANIMATED_PROJECTION`. */
  float wander = 0.0f;    /**< Fraction of the projection random-walk delta
                                absorbed each frame. */
  float camera_wander = 0.0f;    /**< Same, for the outer camera random walk. */
  float central_meridian = 0.0f; /**< Central meridian handed to projections
                                      that take one, in radians. */

  /** @brief Per-parameter id, member, label, range, curve and gate. */
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
          "Central Meridian", 0.0f, math::TWO_PI_F,
          FieldCurve::SHORTEST_PERIODIC, FieldGate::CENTRAL_MERIDIAN},
  };
};
static_assert(field_ids_unique<ProjectionParams>());
static_assert(field_defaults_in_range<ProjectionParams>());

/** @brief Half of the sphere a `Gnomonic` projection images; FOLDED
    overlays both halves. */
enum class GnomonicHemisphere : uint8_t { FOLDED, FRONT, BACK };

/** @brief `ProjectionProvenance::flags` bit set by folded projections. */
inline constexpr uint8_t FOLDED_FLAG = 1U << 0;
/** Render-space divisor floor that caps gnomonic coordinates near 1000. */
inline constexpr float GNOMONIC_AXIS_EPS = 1e-3f;

/**
 * @brief Rational singularity weight f^2 r / (f^2 r + s), with r and s the
 *        squared distances and f the fade sharpness floored at 1e-3.
 * @param regular_distance_sq Squared-distance term that vanishes at the
 *                            singularity.
 * @param singular_distance_sq Squared-distance term that vanishes at the
 *                             regular locus.
 * @param singularity_fade Attenuation sharpness: 1 reaches the regular locus;
 *                         20 confines the fade to a narrow cap.
 * @return Weight in [0, 1]: 0 at the singularity, 1 on the regular locus.
 */
__attribute__((always_inline)) inline float
singularity_attenuation(float regular_distance_sq, float singular_distance_sq,
                        float singularity_fade) {
  const float pf = singularity_fade > 1e-3f ? singularity_fade : 1e-3f;
  const float scaled_distance_sq = pf * pf * regular_distance_sq;
  return scaled_distance_sq / (scaled_distance_sq + singular_distance_sq);
}

/**
 * @brief Singularity weight of the equirectangular projection, singular at the
 *        y = +-1 poles.
 * @param input Unit direction in the projection frame.
 * @param singularity_fade Attenuation sharpness, as in
 *        `singularity_attenuation`.
 * @return Weight in [0, 1]; 0 at the poles.
 */
__attribute__((always_inline)) inline float
equirectangular_weight(const math::Vector &input, float singularity_fade) {
  return singularity_attenuation(input.x * input.x + input.z * input.z,
                                 input.y * input.y, singularity_fade);
}

/**
 * @brief Singularity weight of the Peirce projection, singular at the four
 *        equatorial points midway between the meridian-rotated x and z axes.
 * @param input Unit direction in the projection frame.
 * @param meridian_cos Cosine of the central meridian.
 * @param meridian_sin Sine of the central meridian.
 * @param singularity_fade Attenuation sharpness, as in
 *        `singularity_attenuation`.
 * @param folded Also fade toward the |x| = |z| fold diagonals of the y < 0
 *        hemisphere (DIAMOND and SQUARE layouts).
 * @return Weight in [0, 1]; 0 at a singularity.
 */
__attribute__((always_inline)) inline float
peirce_weight(const math::Vector &input, float meridian_cos, float meridian_sin,
              float singularity_fade, bool folded) {
  const float rotated_x = input.x * meridian_cos + input.z * meridian_sin;
  const float rotated_z = input.z * meridian_cos - input.x * meridian_sin;
  const float singular_cosine =
      (fabsf(rotated_z) + fabsf(rotated_x)) * 0.7071067811865475f;
  float sin_distance_sq = fmaxf(0.0f, 1.0f - singular_cosine * singular_cosine);
  if (folded && input.y < 0.0f) {
    const float fold_sine =
        fabsf(fabsf(rotated_z) - fabsf(rotated_x)) * 0.7071067811865475f;
    sin_distance_sq = fminf(sin_distance_sq, fold_sine * fold_sine);
  }
  return singularity_attenuation(
      sin_distance_sq, fmaxf(0.0f, 1.0f - sin_distance_sq), singularity_fade);
}

/** @brief `peirce_weight` for the folded (DIAMOND/SQUARE) layouts.
    @param input Unit direction in the projection frame.
    @param central_meridian Central meridian, in radians.
    @param singularity_fade Attenuation sharpness, as in
           `singularity_attenuation`.
    @return Weight in [0, 1]; 0 at a singularity. */
__attribute__((always_inline)) inline float
peirce_folded_weight(const math::Vector &input, float central_meridian,
                     float singularity_fade) {
  return peirce_weight(input, cosf(central_meridian), sinf(central_meridian),
                       singularity_fade, true);
}

/**
 * @brief Stereographic projection, singular at the +Y pole.
 * @param input Unit direction in the projection frame.
 * @param singularity_fade Attenuation sharpness, as in
 *        `singularity_attenuation`.
 * @return Plane coordinates and provenance; edge distance is 1 - y.
 */
__attribute__((always_inline)) inline ProjectionResult
stereographic(const math::Vector &input, float singularity_fade) {
  const math::Complex coords = projections::stereo(input);
  return {coords,
          {.region_id = 0,
           .component_id = 0,
           .boundary_flags = static_cast<uint8_t>(ProjectionBoundary::SINGULAR),
           .fade_edge_distance = fmaxf(0.0f, 1.0f - input.y),
           .value_weight = singularity_attenuation(fmaxf(0.0f, 1.0f - input.y),
                                                   fmaxf(0.0f, 1.0f + input.y),
                                                   singularity_fade),
           .flags = 0,
           .traits = projections::projection_traits(
               projections::ProjectionTrait::SINGULAR)}};
}

/**
 * @brief Sinusoidal projection folded about the central meridian.
 * @param input Unit direction in the projection frame.
 * @param central_meridian Fold meridian, in radians.
 * @return Plane coordinates and provenance; `region_id` is 1 on the
 *         negative-longitude half. No edge distance.
 */
#if defined(__EMSCRIPTEN__)
__attribute__((noinline))
#else
__attribute__((always_inline))
#endif
inline ProjectionResult folded_sinusoidal(const math::Vector &input,
                                          float central_meridian) {
  float longitude;
  const math::Complex coords =
      projections::folded_sinusoidal(input, central_meridian, &longitude);
  return {coords,
          {.region_id = static_cast<uint8_t>(longitude < 0.0f),
           .component_id = 0,
           .boundary_flags = 0,
           .fade_edge_distance = projections::NO_EDGE_DISTANCE,
           .value_weight = 1.0f,
           .flags = FOLDED_FLAG,
           .traits = projections::projection_traits(
               projections::ProjectionTrait::FOLDED)}};
}

/**
 * @brief Equirectangular projection cut at the antimeridian.
 * @param input Unit direction in the projection frame.
 * @param central_meridian Longitude at the image centre, in radians.
 * @param singularity_fade Attenuation sharpness, as in
 *        `singularity_attenuation`.
 * @return Plane coordinates and provenance; edge distance is pi - |re|.
 */
__attribute__((always_inline)) inline ProjectionResult
equirectangular(const math::Vector &input, float central_meridian,
                float singularity_fade) {
  const math::Complex coords =
      projections::equirectangular(input, central_meridian);
  return {coords,
          {.region_id = 0,
           .component_id = 0,
           .boundary_flags = static_cast<uint8_t>(ProjectionBoundary::CUT),
           .fade_edge_distance = math::PI_F - fabsf(coords.re),
           .value_weight = equirectangular_weight(input, singularity_fade),
           .flags = 0,
           .traits = projections::projection_traits(
               projections::ProjectionTrait::CUT,
               projections::ProjectionTrait::SINGULAR)}};
}

/**
 * @brief Gnomonic projection (x / y, z / y), singular on the y = 0 circle.
 * @param input Unit direction in the projection frame; |y| is floored at
 *        `GNOMONIC_AXIS_EPS`.
 * @param singularity_fade Attenuation sharpness, as in
 *        `singularity_attenuation`.
 * @param hemisphere Half kept; `domain_coverage` is 0 outside it.
 * @return Plane coordinates and provenance; `region_id` and `component_id`
 *         are 1 for y < 0.
 */
__attribute__((always_inline)) inline ProjectionResult
gnomonic(const math::Vector &input, float singularity_fade,
         GnomonicHemisphere hemisphere) {
  float y = input.y;
  if (fabsf(y) < GNOMONIC_AXIS_EPS)
    y = y < 0.0f ? -GNOMONIC_AXIS_EPS : GNOMONIC_AXIS_EPS;
  const math::Complex coords(input.x / y, input.z / y);
  const bool in_domain =
      hemisphere == GnomonicHemisphere::FOLDED ||
      (hemisphere == GnomonicHemisphere::FRONT ? input.y >= 0.0f
                                               : input.y < 0.0f);
  return {coords,
          {.region_id = static_cast<uint8_t>(input.y < 0.0f),
           .component_id = static_cast<uint8_t>(input.y < 0.0f),
           .boundary_flags = static_cast<uint8_t>(
               static_cast<uint8_t>(ProjectionBoundary::CUT) |
               static_cast<uint8_t>(ProjectionBoundary::SINGULAR)),
           .fade_edge_distance = fabsf(input.y),
           .value_weight = singularity_attenuation(
               input.y * input.y, input.x * input.x + input.z * input.z,
               singularity_fade),
           .flags = 0,
           .traits = projections::projection_traits(
               projections::ProjectionTrait::CUT,
               projections::ProjectionTrait::SINGULAR,
               hemisphere == GnomonicHemisphere::FOLDED
                   ? projections::ProjectionTrait::FOLDED
                   : projections::ProjectionTrait::NONE),
           .edge_class = 0,
           .domain_coverage = in_domain ? 1.0f : 0.0f}};
}

/**
 * @brief Wraps a kernel result as a ProjectionResult.
 * @param result Kernel output.
 * @param coordinate_scale Factor applied to the coordinates; its magnitude
 *        scales `fade_edge_distance`.
 * @param value_weight Weight stored in the provenance.
 * @return The scaled result with the kernel's provenance.
 */
__attribute__((always_inline)) inline ProjectionResult
from_kernel(const projections::ProjectionKernelResult &result,
            float coordinate_scale, float value_weight = 1.0f) {
  return {{result.coords.re * coordinate_scale,
           result.coords.im * coordinate_scale},
          {.region_id = result.region_id,
           .component_id = result.component_id,
           .boundary_flags = result.boundary_flags,
           .fade_edge_distance =
               result.fade_edge_distance * fabsf(coordinate_scale),
           .value_weight = value_weight,
           .flags = result.flags,
           .traits = result.traits,
           .edge_class = result.edge_class}};
}

/**
 * @brief Bonne pseudoconical equal-area projection.
 * @param input Unit direction in the projection frame.
 * @param central_meridian Longitude of the image axis, in radians.
 * @param standard_parallel Signed standard parallel, in radians.
 * @param coordinate_scale Factor applied to the plane coordinates.
 * @return Plane coordinates and provenance.
 */
__attribute__((always_inline)) inline ProjectionResult
bonne(const math::Vector &input, float central_meridian,
      float standard_parallel, float coordinate_scale) {
  return from_kernel(
      projections::bonne_projection(input, central_meridian, standard_parallel),
      coordinate_scale);
}

/**
 * @brief Peirce quincuncial projection with a precomputed meridian rotation.
 * @param input Unit direction in the projection frame.
 * @param central_meridian Longitude of the image axis, in radians.
 * @param layout A `projections::PeirceLayout` value.
 * @param layout_scroll Fraction of a period to translate a strip layout by.
 * @param edge_distance_required Measure `fade_edge_distance`; when false it
 *        is `projections::NO_EDGE_DISTANCE`.
 * @param coordinate_scale Factor applied to the plane coordinates.
 * @param singularity_fade Attenuation sharpness, as in
 *        `singularity_attenuation`.
 * @param meridian_cos Cosine of @p central_meridian.
 * @param meridian_sin Sine of @p central_meridian.
 * @return Plane coordinates and provenance.
 */
__attribute__((always_inline)) inline ProjectionResult
peirce(const math::Vector &input, float central_meridian, uint8_t layout,
       float layout_scroll, bool edge_distance_required, float coordinate_scale,
       float singularity_fade, float meridian_cos, float meridian_sin) {
  return from_kernel(
      projections::peirce_projection(
          input, central_meridian,
          static_cast<projections::PeirceLayout>(layout), layout_scroll,
          edge_distance_required),
      coordinate_scale,
      peirce_weight(
          input, meridian_cos, meridian_sin, singularity_fade,
          layout == static_cast<uint8_t>(projections::PeirceLayout::DIAMOND) ||
              layout ==
                  static_cast<uint8_t>(projections::PeirceLayout::SQUARE)));
}

/**
 * @brief Peirce quincuncial projection.
 * @param input Unit direction in the projection frame.
 * @param central_meridian Longitude of the image axis, in radians.
 * @param layout A `projections::PeirceLayout` value.
 * @param layout_scroll Fraction of a period to translate a strip layout by.
 * @param edge_distance_required Measure `fade_edge_distance`; when false it
 *        is `projections::NO_EDGE_DISTANCE`.
 * @param coordinate_scale Factor applied to the plane coordinates.
 * @param singularity_fade Attenuation sharpness, as in
 *        `singularity_attenuation`.
 * @return Plane coordinates and provenance.
 */
__attribute__((always_inline)) inline ProjectionResult
peirce(const math::Vector &input, float central_meridian, uint8_t layout,
       float layout_scroll, bool edge_distance_required, float coordinate_scale,
       float singularity_fade) {
  return peirce(input, central_meridian, layout, layout_scroll,
                edge_distance_required, coordinate_scale, singularity_fade,
                cosf(central_meridian), sinf(central_meridian));
}

/**
 * @brief Approximate square-layout Peirce projection at a zero central
 *        meridian.
 * @param input Unit direction in the projection frame.
 * @param coordinate_scale Factor applied to the plane coordinates.
 * @param singularity_fade Attenuation sharpness, as in
 *        `singularity_attenuation`.
 * @return Plane coordinates and provenance.
 */
__attribute__((always_inline)) inline ProjectionResult
peirce_fast_square(const math::Vector &input, float coordinate_scale,
                   float singularity_fade) {
  return from_kernel(projections::peirce_projection_fast_square(input),
                     coordinate_scale,
                     peirce_folded_weight(input, 0.0f, singularity_fade));
}

template <typename State, typename Binding>
concept FrameProvider = Detail::ProviderFor<State, Binding> &&
                        requires(const typename Binding::FrameState &frame) {
                          {
                            State::conjugate(frame)
                          } -> std::same_as<const math::Quaternion &>;
                        };

/** @brief Exact projection base whose frame conjugate is the provider
    @p State's. */
template <typename State>
struct FrameConjugateFromState : ApproximationDefaults {
  /**
   * @brief Rotation taking world directions into the projection frame.
   * @param frame Frame state.
   * @return `State`'s frame conjugate.
   */
  __attribute__((always_inline)) static const math::Quaternion &
  frame_conjugate(const typename State::FrameState &frame) {
    return State::conjugate(frame);
  }
};

/**
 * @brief Airocean icosahedral net.
 * @param input Unit direction in the projection frame.
 * @param central_meridian Longitude of the net's axis, in radians.
 * @param horizontal Turn the finished net a quarter turn.
 * @param edge_distance_required Measure the per-edge cut distances.
 * @param coordinate_scale Factor applied to the plane coordinates.
 * @return Plane coordinates and provenance.
 */
__attribute__((always_inline)) inline ProjectionResult
airocean(const math::Vector &input, float central_meridian, bool horizontal,
         bool edge_distance_required, float coordinate_scale) {
  return from_kernel(
      projections::airocean_projection_meridian(
          input, central_meridian, horizontal, edge_distance_required),
      coordinate_scale);
}

/**
 * @brief Airocean icosahedral net with a precomputed meridian rotation.
 * @param input Unit direction in the projection frame.
 * @param horizontal Turn the finished net a quarter turn.
 * @param edge_distance_required Measure the per-edge cut distances.
 * @param coordinate_scale Factor applied to the plane coordinates.
 * @param meridian_cos Cosine of the central meridian.
 * @param meridian_sin Sine of the central meridian.
 * @return Plane coordinates and provenance.
 */
__attribute__((always_inline)) inline ProjectionResult
airocean(const math::Vector &input, bool horizontal,
         bool edge_distance_required, float coordinate_scale,
         float meridian_cos, float meridian_sin) {
  return from_kernel(projections::airocean_projection(input, meridian_cos,
                                                      meridian_sin, horizontal,
                                                      edge_distance_required),
                     coordinate_scale);
}

/** @brief Bonne pseudoconical equal-area projection; `North` picks the sign of
    the standard parallel, and so the hemisphere the cone opens toward. */
template <typename State, bool North>
struct Bonne : FrameConjugateFromState<State> {
  /// `fade_edge_distance` is always measured.
  static constexpr bool EDGE_DISTANCE_AVAILABLE = true;
  /// Binding of the frame-state provider `State`.
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the frame
   *        conjugate, central meridian, standard parallel and
   *        coordinate scale.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::standard_parallel(frame) } -> std::same_as<float>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
      };

  /**
   * @brief Projects a direction given in the projection frame.
   * @param input Unit direction in the projection frame.
   * @param frame Frame state.
   * @return Plane coordinates and provenance.
   */
  __attribute__((always_inline)) static ProjectionResult
  project(const math::Vector &input, const FrameState &frame) {
    const float hemisphere = North ? 1.0f : -1.0f;
    return bonne(input, State::central_meridian(frame),
                 hemisphere * State::standard_parallel(frame),
                 State::coordinate_scale(frame));
  }
};

/** @brief Stereographic projection: conformal, with one singular pole the
    singularity fade attenuates. */
template <typename State>
struct Stereographic : FrameConjugateFromState<State> {
  /// `fade_edge_distance` is always measured.
  static constexpr bool EDGE_DISTANCE_AVAILABLE = true;
  /// Binding of the frame-state provider `State`.
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the frame
   *        conjugate and singularity fade.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::singularity_fade(frame) } -> std::same_as<float>;
      };

  /**
   * @brief Projects a direction given in the projection frame.
   * @param input Unit direction in the projection frame.
   * @param frame Frame state.
   * @return Plane coordinates and provenance.
   */
  __attribute__((always_inline)) static ProjectionResult
  project(const math::Vector &input, const FrameState &frame) {
    return stereographic(input, State::singularity_fade(frame));
  }
};

/** @brief Sinusoidal projection with the azimuth folded about the central
    meridian: both hemispheres share one image, and there is no singular locus
    to attenuate. */
template <typename State>
struct FoldedSinusoidal : FrameConjugateFromState<State> {
  /// Binding of the frame-state provider `State`.
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the frame
   *        conjugate and central meridian.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
      };

  /**
   * @brief Projects a direction given in the projection frame.
   * @param input Unit direction in the projection frame.
   * @param frame Frame state.
   * @return Plane coordinates and provenance.
   */
  __attribute__((always_inline)) static ProjectionResult
  project(const math::Vector &input, const FrameState &frame) {
    return folded_sinusoidal(input, State::central_meridian(frame));
  }
};

/** @brief Equirectangular projection: cut at the antimeridian, with both
    poles attenuated by the singularity fade. */
template <typename State>
struct Equirectangular : FrameConjugateFromState<State> {
  /// `fade_edge_distance` is always measured.
  static constexpr bool EDGE_DISTANCE_AVAILABLE = true;
  /// Binding of the frame-state provider `State`.
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the frame
   *        conjugate, central meridian and singularity fade.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::singularity_fade(frame) } -> std::same_as<float>;
      };

  /**
   * @brief Projects a direction given in the projection frame.
   * @param input Unit direction in the projection frame.
   * @param frame Frame state.
   * @return Plane coordinates and provenance.
   */
  __attribute__((always_inline)) static ProjectionResult
  project(const math::Vector &input, const FrameState &frame) {
    return equirectangular(input, State::central_meridian(frame),
                           State::singularity_fade(frame));
  }
};

/** @brief Gnomonic projection about the Y axis, singular on the y = 0 great
    circle; `Hemisphere` folds the two halves together or keeps one. */
template <typename State, GnomonicHemisphere Hemisphere>
struct Gnomonic : FrameConjugateFromState<State> {
  /// `fade_edge_distance` is always measured.
  static constexpr bool EDGE_DISTANCE_AVAILABLE = true;
  /// Binding of the frame-state provider `State`.
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the frame
   *        conjugate and singularity fade.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::singularity_fade(frame) } -> std::same_as<float>;
      };

  /**
   * @brief Projects a direction given in the projection frame.
   * @param input Unit direction in the projection frame.
   * @param frame Frame state.
   * @return Plane coordinates and provenance.
   */
  __attribute__((always_inline)) static ProjectionResult
  project(const math::Vector &input, const FrameState &frame) {
    return gnomonic(input, State::singularity_fade(frame), Hemisphere);
  }
};

/** @brief Central-meridian rotation prepared once per frame. */
struct PreparedMeridian {
  float cosine; ///< Cosine of the central meridian.
  float sine;   ///< Sine of the central meridian.

  /**
   * @brief Prepares the rotation for one central meridian.
   * @param angle Central meridian, in radians.
   * @return Cosine and sine of @p angle.
   */
  static PreparedMeridian from_angle(float angle) {
    return {cosf(angle), sinf(angle)};
  }
};

/** @brief Peirce quincuncial projection, conformal but for four singularities;
    `Layout` picks diamond, square or strip tiling and `EdgeDistanceRequired`
    makes the kernel compute edge distance unconditionally. */
template <typename State, uint8_t Layout, bool EdgeDistanceRequired>
struct Peirce : FrameConjugateFromState<State> {
  /// Binding of the frame-state provider `State`.
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.
  using Prepared = PreparedMeridian; ///< Per-frame central-meridian rotation.

  /**
   * @brief Caches the central meridian's cosine and sine for the frame.
   * @param frame Frame state.
   * @return The prepared meridian rotation.
   */
  static Prepared prepare(const FrameState &frame) {
    return Prepared::from_angle(State::central_meridian(frame));
  }

  /// Whether `fade_edge_distance` is measured.
  static constexpr bool EDGE_DISTANCE_AVAILABLE = EdgeDistanceRequired;

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the frame
   *        conjugate, central meridian, layout scroll, coordinate
   *        scale and singularity fade.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::layout_scroll(frame) } -> std::same_as<float>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
        { State::singularity_fade(frame) } -> std::same_as<float>;
      };

  /**
   * @brief Projects a direction given in the projection frame.
   * @param input Unit direction in the projection frame.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Plane coordinates and provenance.
   */
  __attribute__((always_inline)) static ProjectionResult
  project(const math::Vector &input, const FrameState &frame,
          const Prepared &prepared) {
    return peirce(input, State::central_meridian(frame), Layout,
                  State::layout_scroll(frame), EdgeDistanceRequired,
                  State::coordinate_scale(frame),
                  State::singularity_fade(frame), prepared.cosine,
                  prepared.sine);
  }
  /**
   * @brief Projects a direction, preparing the meridian rotation inline.
   * @param input Unit direction in the projection frame.
   * @param frame Frame state.
   * @return Plane coordinates and provenance.
   */
  static ProjectionResult project(const math::Vector &input,
                                  const FrameState &frame) {
    return project(input, frame, prepare(frame));
  }
};

/** @brief Approximation bounds of the fast square Peirce path. */
inline constexpr std::array<ApproximationMetric, 3> PEIRCE_FAST_SQUARE_METRICS{{
    {ApproximationDomain::PROJECTED_COORDINATE,
     ApproximationAggregation::MAXIMUM, 1.2e-3f, "plane units"},
    {ApproximationDomain::PROJECTED_EDGE_DISTANCE,
     ApproximationAggregation::MAXIMUM, 2e-4f,
     "radians * abs(coordinate_scale)"},
    {ApproximationDomain::FRAMEBUFFER, ApproximationAggregation::MAXIMUM,
     256.0f, "channel code"},
}};

/** @brief Approximate square-layout Peirce projection; the provider must pin
    the central meridian to zero. */
template <typename State>
struct PeirceFastSquare : FrameConjugateFromState<State> {
  /// Binding of the frame-state provider `State`.
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /// `fade_edge_distance` is always measured.
  static constexpr bool EDGE_DISTANCE_AVAILABLE = true;
  /// Takes the approximate fast-square kernel.
  static constexpr bool APPROXIMATE = true;
  /// Reference the approximation is checked against.
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::PEIRCE_FAST_SQUARE;
  static constexpr std::array<ApproximationMetric, 3> METRICS =
      PEIRCE_FAST_SQUARE_METRICS; ///< Approximation error bounds.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the frame
   *        conjugate, coordinate scale and singularity fade, and
   *        pins the central meridian to zero.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::ZERO_CENTRAL_MERIDIAN } -> std::convertible_to<bool>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
        { State::singularity_fade(frame) } -> std::same_as<float>;
      } && State::ZERO_CENTRAL_MERIDIAN;

  /**
   * @brief Projects a direction given in the projection frame.
   * @param input Unit direction in the projection frame.
   * @param frame Frame state.
   * @return Plane coordinates and provenance.
   */
  __attribute__((always_inline)) static ProjectionResult
  project(const math::Vector &input, const FrameState &frame) {
    return peirce_fast_square(input, State::coordinate_scale(frame),
                              State::singularity_fade(frame));
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
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.
  using Prepared = PreparedMeridian; ///< Per-frame central-meridian rotation.

  /**
   * @brief Caches the central meridian's cosine and sine for the frame.
   * @param frame Frame state.
   * @return The prepared meridian rotation.
   */
  static Prepared prepare(const FrameState &frame) {
    return Prepared::from_angle(State::central_meridian(frame));
  }

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the frame
   *        conjugate, central meridian, layout scroll, coordinate
   *        scale and singularity fade.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::layout_scroll(frame) } -> std::same_as<float>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
        { State::singularity_fade(frame) } -> std::same_as<float>;
      };

  /**
   * @brief Projects a direction given in the projection frame.
   * @param input Unit direction in the projection frame.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Plane coordinates and provenance.
   */
  __attribute__((always_inline)) static ProjectionResult
  project(const math::Vector &input, const FrameState &frame,
          const Prepared &prepared) {
    if (State::central_meridian(frame) == 0.0f)
      return peirce_fast_square(input, State::coordinate_scale(frame),
                                State::singularity_fade(frame));
    return peirce(
        input, State::central_meridian(frame),
        static_cast<uint8_t>(projections::PeirceLayout::SQUARE),
        State::layout_scroll(frame), true, State::coordinate_scale(frame),
        State::singularity_fade(frame), prepared.cosine, prepared.sine);
  }
  /**
   * @brief Projects a direction, preparing the meridian rotation inline.
   * @param input Unit direction in the projection frame.
   * @param frame Frame state.
   * @return Plane coordinates and provenance.
   */
  static ProjectionResult project(const math::Vector &input,
                                  const FrameState &frame) {
    return project(input, frame, prepare(frame));
  }
};

/** @brief Airocean icosahedral net; `Horizontal` turns the finished net a
    quarter turn and `EdgeDistanceRequired` makes the kernel compute the
    per-edge cut distances unconditionally. */
template <typename State, bool Horizontal, bool EdgeDistanceRequired>
struct Airocean : FrameConjugateFromState<State> {
  /// Binding of the frame-state provider `State`.
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.
  using Prepared = PreparedMeridian; ///< Per-frame central-meridian rotation.

  /**
   * @brief Caches the central meridian's cosine and sine for the frame.
   * @param frame Frame state.
   * @return The prepared meridian rotation.
   */
  static Prepared prepare(const FrameState &frame) {
    return Prepared::from_angle(State::central_meridian(frame));
  }

  /// Whether `fade_edge_distance` is measured.
  static constexpr bool EDGE_DISTANCE_AVAILABLE = EdgeDistanceRequired;

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the frame
   *        conjugate, central meridian and coordinate scale.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
      };

  /**
   * @brief Projects a direction given in the projection frame.
   * @param input Unit direction in the projection frame.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Plane coordinates and provenance.
   */
  __attribute__((always_inline)) static ProjectionResult
  project(const math::Vector &input, const FrameState &frame,
          const Prepared &prepared) {
    return airocean(input, Horizontal, EdgeDistanceRequired,
                    State::coordinate_scale(frame), prepared.cosine,
                    prepared.sine);
  }
  /**
   * @brief Projects a direction, preparing the meridian rotation inline.
   * @param input Unit direction in the projection frame.
   * @param frame Frame state.
   * @return Plane coordinates and provenance.
   */
  static ProjectionResult project(const math::Vector &input,
                                  const FrameState &frame) {
    return project(input, frame, prepare(frame));
  }
};

} // namespace Projection

} // namespace Pullback
