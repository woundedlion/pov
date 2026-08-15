/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <tuple>
#include <type_traits>
#include <utility>

#include "color/color.h"
#include "engine/memory.h"
#include "math/geometry.h"
#include "math/lenses.h"
#include "math/noise_field.h"
#include "math/projections.h"
#include "math/stereographic.h"

/**
 * @file pullback.h
 * @brief Typed pullback composition and standard carrier vocabulary.
 */

namespace Pullback {

enum class StageKind : uint8_t {
  OUTER_CAMERA,
  SURFACE_PROJECT,
  PLANAR_WARP,
  SOURCE,
  MATERIAL,
  COLOR
};

enum class CodeEmission : uint8_t {
  INLINE_ONLY,
  OUT_OF_LINE_FLASH,
  OUT_OF_LINE_ITCM
};

enum class ApproximationOracleId : uint8_t {
  NONE,
  PEIRCE_FAST_SQUARE,
  HUE_ROTATION_AND_NOISE_LUTS
};

enum class ApproximationDomain : uint8_t {
  PROJECTED_COORDINATE,
  PROJECTED_EDGE_DISTANCE,
  COLOR_CHANNEL,
  FRAMEBUFFER
};

enum class ApproximationAggregation : uint8_t { MAXIMUM, MEAN };

struct ApproximationMetric {
  ApproximationDomain domain;
  ApproximationAggregation aggregation;
  float limit;
  const char *unit;
};

enum class ProjectionBoundary : uint8_t { CUT = 1U << 0, SINGULAR = 1U << 1 };

struct ProjectionSample {
  Complex coords;
  uint8_t region_id;
  uint8_t component_id;
  uint8_t boundary_flags;
  float fade_edge_distance;
  float value_weight;
  uint8_t flags;
  uint8_t traits = 0;
  uint8_t edge_class = 0;
  float domain_coverage = 1.0f;
  Vector sphere = Vector();
  float surface_path_length = 0.0f;
};

struct SurfaceResult {
  Vector sphere;
  float path_length;
};

struct WarpStepResult {
  Complex coords;
  Complex delta;
  float deformation;
  float path_length;
};

struct WarpResult {
  Complex coords;
  Complex net_delta;
  float path_length;
};

struct SourceInput {
  ProjectionSample projected;
  WarpResult warped;
};

struct MaterialInput {
  ProjectionSample projected;
  WarpResult warped;
  float field;
};

struct MaterialSample {
  float value;
  float coverage;
  Vector sphere;
  float path_length;
};

enum class ProfileEvent : uint8_t {
  LENS,
  SURFACE_NOISE,
  PROJECTION,
  PLANAR_WARP,
  MIRROR_TILE,
  SOURCE,
  MATERIAL,
  COLOR
};

struct NoInstrumentation {
  struct Token {};

  __attribute__((always_inline)) static Token mark() { return {}; }

  template <ProfileEvent Event>
  __attribute__((always_inline)) static void span(Token) {
    static_cast<void>(Event);
  }
};

template <typename Stage, typename Binding, typename = void>
struct HasStageContract : std::false_type {};

template <typename Stage, typename Binding>
struct HasStageContract<
    Stage, Binding,
    std::void_t<typename Stage::Binding, typename Stage::FrameState,
                typename Stage::Input, typename Stage::Output,
                decltype(Stage::KIND), decltype(Stage::EMISSION),
                decltype(Stage::APPROXIMATE), decltype(Stage::TERMINAL),
                decltype(Stage::ORACLE), decltype(Stage::METRICS),
                decltype(Stage::NON_FLOATING_FIELDS_EXACT),
                decltype(Stage::run(
                    std::declval<const typename Stage::Input &>(),
                    std::declval<const typename Stage::FrameState &>()))>>
    : std::true_type {};

template <typename Stage, typename Binding,
          bool Present = HasStageContract<Stage, Binding>::value>
struct HasTypedStageContract : std::false_type {};

template <typename T> struct IsApproximationMetricArray : std::false_type {};

template <size_t Size>
struct IsApproximationMetricArray<std::array<ApproximationMetric, Size>>
    : std::true_type {};

template <typename Stage, typename Binding>
struct HasTypedStageContract<Stage, Binding, true>
    : std::bool_constant<
          std::is_same_v<typename Stage::FrameState,
                         typename Stage::Binding::FrameState> &&
          std::is_same_v<decltype(Stage::KIND), const StageKind> &&
          std::is_same_v<decltype(Stage::EMISSION), const CodeEmission> &&
          std::is_same_v<decltype(Stage::APPROXIMATE), const bool> &&
          std::is_same_v<decltype(Stage::TERMINAL), const bool> &&
          std::is_same_v<decltype(Stage::NON_FLOATING_FIELDS_EXACT),
                         const bool> &&
          std::is_same_v<decltype(Stage::ORACLE),
                         const ApproximationOracleId> &&
          IsApproximationMetricArray<
              std::remove_cv_t<decltype(Stage::METRICS)>>::value> {};

template <typename Stage> consteval bool has_final_framebuffer_metric() {
  for (const ApproximationMetric &metric : Stage::METRICS)
    if (metric.domain == ApproximationDomain::FRAMEBUFFER)
      return true;
  return false;
}

template <typename Binding, typename... Stages>
consteval bool binding_extra_validation() {
  if constexpr (requires {
                  Binding::template ExtraValidation<Stages...>::value;
                })
    return Binding::template ExtraValidation<Stages...>::value;
  else
    return true;
}

template <bool LevelOneValid, typename Binding, typename... Stages>
struct PipelineValidationLevel2 {
  static constexpr bool BINDINGS = false;
  static constexpr bool EMPTY_POLICIES = false;
  static constexpr bool ORDER = false;
  static constexpr bool RUN_RETURNS = false;
  static constexpr bool CARRIERS = false;
  static constexpr bool TERMINALS = false;
  static constexpr bool APPROXIMATIONS = false;
  static constexpr bool EXTRA_VALIDATION = false;
};

template <typename Binding, typename... Stages>
struct PipelineValidationLevel2<true, Binding, Stages...> {
  using OuterStage = std::tuple_element_t<0, std::tuple<Stages...>>;
  using SurfaceStage = std::tuple_element_t<1, std::tuple<Stages...>>;
  using WarpStage = std::tuple_element_t<2, std::tuple<Stages...>>;
  using SourceStage = std::tuple_element_t<3, std::tuple<Stages...>>;
  using MaterialStage = std::tuple_element_t<4, std::tuple<Stages...>>;
  using ColorStage = std::tuple_element_t<5, std::tuple<Stages...>>;

  static constexpr bool BINDINGS =
      (std::is_same_v<typename Stages::Binding, Binding> && ...);
  static constexpr bool EMPTY_POLICIES = (std::is_empty_v<Stages> && ...);
  static constexpr bool ORDER =
      OuterStage::KIND == StageKind::OUTER_CAMERA &&
      SurfaceStage::KIND == StageKind::SURFACE_PROJECT &&
      WarpStage::KIND == StageKind::PLANAR_WARP &&
      SourceStage::KIND == StageKind::SOURCE &&
      MaterialStage::KIND == StageKind::MATERIAL &&
      ColorStage::KIND == StageKind::COLOR;
  static constexpr bool RUN_RETURNS =
      (std::is_same_v<decltype(Stages::run(
                          std::declval<const typename Stages::Input &>(),
                          std::declval<const typename Stages::FrameState &>())),
                      typename Stages::Output> &&
       ...);
  static constexpr bool CARRIERS =
      std::is_same_v<typename OuterStage::Input, Vector> &&
      std::is_same_v<typename OuterStage::Output,
                     typename SurfaceStage::Input> &&
      std::is_same_v<typename SurfaceStage::Output,
                     typename WarpStage::Input> &&
      std::is_same_v<typename WarpStage::Output, typename SourceStage::Input> &&
      std::is_same_v<typename SourceStage::Output,
                     typename MaterialStage::Input> &&
      std::is_same_v<typename MaterialStage::Output,
                     typename ColorStage::Input> &&
      std::is_same_v<typename ColorStage::Output, Color4>;
  static constexpr bool TERMINALS =
      !OuterStage::TERMINAL && !SurfaceStage::TERMINAL &&
      !WarpStage::TERMINAL && !SourceStage::TERMINAL &&
      !MaterialStage::TERMINAL && ColorStage::TERMINAL;
  static constexpr bool APPROXIMATIONS =
      (((Stages::APPROXIMATE && Stages::ORACLE != ApproximationOracleId::NONE &&
         Stages::METRICS.size() != 0 && Stages::NON_FLOATING_FIELDS_EXACT &&
         has_final_framebuffer_metric<Stages>()) ||
        (!Stages::APPROXIMATE &&
         Stages::ORACLE == ApproximationOracleId::NONE &&
         Stages::METRICS.size() == 0)) &&
       ...);
  static constexpr bool EXTRA_VALIDATION =
      binding_extra_validation<Binding, Stages...>();
};

template <typename Binding, typename... Stages>
struct PipelineValidation
    : PipelineValidationLevel2<
          sizeof...(Stages) == 6 &&
              (HasTypedStageContract<Stages, Binding>::value && ...),
          Binding, Stages...> {
  static constexpr bool ARITY = sizeof...(Stages) == 6;
  static constexpr bool CONTRACTS =
      (HasTypedStageContract<Stages, Binding>::value && ...);
  using Level2 =
      PipelineValidationLevel2<ARITY && CONTRACTS, Binding, Stages...>;

  static constexpr bool BINDINGS = Level2::BINDINGS;
  static constexpr bool EMPTY_POLICIES = Level2::EMPTY_POLICIES;
  static constexpr bool ORDER = Level2::ORDER;
  static constexpr bool RUN_RETURNS = Level2::RUN_RETURNS;
  static constexpr bool CARRIERS = Level2::CARRIERS;
  static constexpr bool TERMINALS = Level2::TERMINALS;
  static constexpr bool APPROXIMATIONS = Level2::APPROXIMATIONS;
  static constexpr bool EXTRA_VALIDATION = Level2::EXTRA_VALIDATION;
};

template <typename Stage, typename Key>
concept StageMatchesKey = requires(const Key &key) {
  { Stage::implements(key) } -> std::convertible_to<bool>;
};

template <typename BindingT, typename... Stages> struct Pipeline {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  using Validation = PipelineValidation<Binding, Stages...>;

  static constexpr size_t STAGE_COUNT = sizeof...(Stages);
  template <size_t Index>
  using stage_at = std::tuple_element_t<Index, std::tuple<Stages...>>;

  static_assert(Validation::ARITY, "pullback pipeline: wrong stage count");
  static_assert(Validation::CONTRACTS,
                "pullback pipeline: missing or mistyped stage contract");
  static_assert(!Validation::ARITY || !Validation::CONTRACTS ||
                    Validation::BINDINGS,
                "pullback pipeline: stage binding mismatch");
  static_assert(!Validation::ARITY || !Validation::CONTRACTS ||
                    Validation::EMPTY_POLICIES,
                "pullback pipeline: stage policies must be empty");
  static_assert(!Validation::ARITY || !Validation::CONTRACTS ||
                    Validation::ORDER,
                "pullback pipeline: wrong stage order");
  static_assert(!Validation::ARITY || !Validation::CONTRACTS ||
                    Validation::RUN_RETURNS,
                "pullback pipeline: wrong stage return type");
  static_assert(!Validation::ARITY || !Validation::CONTRACTS ||
                    Validation::CARRIERS,
                "pullback pipeline: carrier mismatch");
  static_assert(!Validation::ARITY || !Validation::CONTRACTS ||
                    Validation::TERMINALS,
                "pullback pipeline: terminal contract violation");
  static_assert(!Validation::ARITY || !Validation::CONTRACTS ||
                    Validation::APPROXIMATIONS,
                "pullback pipeline: malformed approximation metadata");
  static_assert(!Validation::ARITY || !Validation::CONTRACTS ||
                    Validation::EXTRA_VALIDATION,
                "pullback pipeline: consumer validation failed");

  __attribute__((always_inline)) static Color4
  evaluate(const Vector &view, const FrameState &frame) {
    return run_stage<0>(view, frame);
  }

  HS_FLASH_MEMBER static Color4 shade(const Vector &view,
                                      const FrameState &frame) {
    return evaluate(view, frame);
  }

  template <typename Key>
    requires(StageMatchesKey<Stages, Key> && ...)
  static constexpr bool implements(const Key &key) {
    return (Stages::implements(key) && ...);
  }

protected:
  template <size_t Index, typename Input>
  __attribute__((always_inline)) static Color4
  run_stage(const Input &input, const FrameState &frame) {
    using CurrentStage = stage_at<Index>;
    static_assert(std::is_same_v<Input, typename CurrentStage::Input>,
                  "pullback pipeline: stage input mismatch");
    const typename CurrentStage::Output output =
        CurrentStage::run(input, frame);
    if constexpr (Index + 1 == STAGE_COUNT)
      return output;
    else
      return run_stage<Index + 1>(output, frame);
  }
};

struct ExactPolicy {
  static constexpr bool APPROXIMATE = false;
  static constexpr bool NON_FLOATING_FIELDS_EXACT = true;
  static constexpr ApproximationOracleId ORACLE = ApproximationOracleId::NONE;
  static constexpr std::array<ApproximationMetric, 0> METRICS{};
};

namespace Detail {

template <typename... Policies> struct FirstApproximate {
  using Type = ExactPolicy;
};

template <typename Head, typename... Tail>
struct FirstApproximate<Head, Tail...> {
  using Type = std::conditional_t<Head::APPROXIMATE, Head,
                                  typename FirstApproximate<Tail...>::Type>;
};

template <typename... Policies> struct CombinedApproximation {
  static constexpr size_t COUNT =
      (static_cast<size_t>(Policies::APPROXIMATE) + ... + 0U);
  static_assert(
      COUNT <= 1,
      "pullback stage: multiple approximation oracles require an explicit "
      "combined oracle");
  using Owner = typename FirstApproximate<Policies...>::Type;
  static constexpr bool APPROXIMATE = Owner::APPROXIMATE;
  static constexpr bool NON_FLOATING_FIELDS_EXACT =
      Owner::NON_FLOATING_FIELDS_EXACT;
  static constexpr ApproximationOracleId ORACLE = Owner::ORACLE;
  static constexpr auto METRICS = Owner::METRICS;
};

template <typename Policy, typename Binding>
consteval bool policy_provider_valid() {
  if constexpr (requires { Policy::template PROVIDER_VALID<Binding>; })
    return Policy::template PROVIDER_VALID<Binding>;
  else
    return true;
}

template <typename BindingT, StageKind KindV, typename InputT, typename OutputT,
          bool TerminalV, typename... Policies>
struct StageContract : CombinedApproximation<Policies...> {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  using Input = InputT;
  using Output = OutputT;

  static constexpr StageKind KIND = KindV;
  static constexpr CodeEmission EMISSION = CodeEmission::INLINE_ONLY;
  static constexpr bool TERMINAL = TerminalV;

  static_assert((std::is_empty_v<Policies> && ...),
                "pullback stage: policies must be empty");
  static_assert((policy_provider_valid<Policies, Binding>() && ...),
                "pullback stage: malformed or foreign provider");
};

template <typename State, typename Binding>
concept ProviderFor =
    std::is_empty_v<State> && std::is_trivially_constructible_v<State> &&
    requires {
      typename State::Binding;
      typename State::FrameState;
    } && std::same_as<typename State::Binding, Binding> &&
    std::same_as<typename State::FrameState, typename Binding::FrameState>;

inline float clamp_unit(float value) {
  if (value <= 0.0f)
    return 0.0f;
  if (value >= 1.0f)
    return 1.0f;
  return value;
}

inline float smooth_ramp(float low, float high, float value) {
  if (low == high)
    return static_cast<float>(value > low);
  const float t = clamp_unit((value - low) / (high - low));
  return t * t * (3.0f - 2.0f * t);
}

} // namespace Detail

namespace Surface {

enum class Integrator : uint8_t { EULER, MIDPOINT, MIDPOINT_2X };

struct Euler {
  static constexpr Integrator VALUE = Integrator::EULER;
};

struct Midpoint {
  static constexpr Integrator VALUE = Integrator::MIDPOINT;
};

struct Midpoint2x {
  static constexpr Integrator VALUE = Integrator::MIDPOINT_2X;
};

struct Identity : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static SurfaceResult
  apply(const Vector &input, const FrameState &) {
    return {input, 0.0f};
  }
};

__attribute__((always_inline)) inline float path_length(const Vector &step,
                                                        bool required) {
  return required ? sqrtf(dot(step, step)) : 0.0f;
}

__attribute__((always_inline)) inline SurfaceResult
finish_step(const Vector &input, const Vector &step,
            bool path_length_required) {
  return {sphere_exp_map_half_radian(input, step),
          path_length(step, path_length_required)};
}

__attribute__((always_inline)) inline SurfaceResult
direct_noise(const Vector &input, const FastNoiseLite &noise,
             ::NoiseBasis basis, float scale, const Vector &loop_offset,
             float strength, float direction_cos, float direction_sin,
             bool path_length_required) {
  if (strength == 0.0f)
    return {input, 0.0f};
  const Vector q = noise_sphere_coordinate(input, scale, loop_offset);
  const Vector tangent = sample_direct_tangent(noise, basis, q, input,
                                               direction_cos, direction_sin);
  return finish_step(input, strength * tangent, path_length_required);
}

__attribute__((always_inline)) inline Vector
curl_field(const Vector &input, const FastNoiseLite &noise, ::NoiseBasis basis,
           float scale, const Vector &loop_offset) {
  const Vector q = noise_sphere_coordinate(input, scale, loop_offset);
  return sample_curl_tangent(noise, basis, q, input);
}

__attribute__((always_inline)) inline SurfaceResult
curl_midpoint_step(const Vector &input, const FastNoiseLite &noise,
                   ::NoiseBasis basis, float scale, const Vector &loop_offset,
                   float distance, bool path_length_required) {
  const Vector first = curl_field(input, noise, basis, scale, loop_offset);
  const Vector midpoint =
      sphere_exp_map_half_radian(input, 0.5f * distance * first);
  const Vector midpoint_field =
      curl_field(midpoint, noise, basis, scale, loop_offset);
  return finish_step(
      input, distance * transport_tangent(midpoint, input, midpoint_field),
      path_length_required);
}

HS_FLASH_INLINE inline SurfaceResult
curl_noise(const Vector &input, const FastNoiseLite &noise, ::NoiseBasis basis,
           Integrator integrator, float scale, const Vector &loop_offset,
           float strength, bool path_length_required) {
  if (strength == 0.0f)
    return {input, 0.0f};
  if (integrator == Integrator::EULER)
    return finish_step(
        input, strength * curl_field(input, noise, basis, scale, loop_offset),
        path_length_required);
  if (integrator == Integrator::MIDPOINT)
    return curl_midpoint_step(input, noise, basis, scale, loop_offset, strength,
                              path_length_required);
  const SurfaceResult first =
      curl_midpoint_step(input, noise, basis, scale, loop_offset,
                         0.5f * strength, path_length_required);
  SurfaceResult second =
      curl_midpoint_step(first.sphere, noise, basis, scale, loop_offset,
                         0.5f * strength, path_length_required);
  second.path_length += first.path_length;
  return second;
}

template <typename State, ::NoiseBasis Basis> struct DirectNoise : ExactPolicy {
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::scale(frame) } -> std::same_as<float>;
        { State::loop_offset(frame) } -> std::same_as<const Vector &>;
        { State::strength(frame) } -> std::same_as<float>;
        { State::direction_cos(frame) } -> std::same_as<float>;
        { State::direction_sin(frame) } -> std::same_as<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      };

  template <typename FrameState>
  __attribute__((always_inline)) static SurfaceResult
  apply(const Vector &input, const FrameState &frame) {
    return direct_noise(input, State::noise(frame), Basis, State::scale(frame),
                        State::loop_offset(frame), State::strength(frame),
                        State::direction_cos(frame),
                        State::direction_sin(frame),
                        State::path_length_required(frame));
  }
};

template <typename State, ::NoiseBasis Basis, typename IntegratorPolicy>
struct CurlNoise : ExactPolicy {
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::scale(frame) } -> std::same_as<float>;
        { State::loop_offset(frame) } -> std::same_as<const Vector &>;
        { State::strength(frame) } -> std::same_as<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      } && requires {
        { IntegratorPolicy::VALUE } -> std::convertible_to<Integrator>;
      };

  template <typename FrameState>
  HS_FLASH_INLINE static SurfaceResult apply(const Vector &input,
                                             const FrameState &frame) {
    return curl_noise(input, State::noise(frame), Basis,
                      IntegratorPolicy::VALUE, State::scale(frame),
                      State::loop_offset(frame), State::strength(frame),
                      State::path_length_required(frame));
  }
};

} // namespace Surface

namespace Lens {

struct Identity : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return input;
  }
};

struct Glitch : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::glitch_lens(input);
  }
};

struct Twist : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::twist_lens(input);
  }
};

struct Kaleidoscope : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::kaleidoscope_lens(input);
  }
};

struct TetrahedralKaleidoscope : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(input,
                                                lenses::TETRAHEDRAL_MIRRORS);
  }
};

struct OctahedralKaleidoscope : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(input,
                                                lenses::OCTAHEDRAL_MIRRORS);
  }
};

struct DodecahedralKaleidoscope : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::dodecahedral_kaleidoscope_lens(input);
  }
};

struct TriangularPrismKaleidoscope : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::TRIANGULAR_PRISM_MIRRORS);
  }
};

struct SquarePrismKaleidoscope : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(input,
                                                lenses::SQUARE_PRISM_MIRRORS);
  }
};

struct PentagonalPrismKaleidoscope : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::PENTAGONAL_PRISM_MIRRORS);
  }
};

struct HexagonalPrismKaleidoscope : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::HEXAGONAL_PRISM_MIRRORS);
  }
};

struct OctagonalPrismKaleidoscope : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::OCTAGONAL_PRISM_MIRRORS);
  }
};

} // namespace Lens

namespace Projection {

enum class GnomonicHemisphere : uint8_t { FOLDED, FRONT, BACK };

static constexpr uint8_t FOLDED_FLAG = 1U << 0;
static constexpr float GNOMONIC_AXIS_EPS = 1e-3f;

__attribute__((always_inline)) inline ProjectionSample
stereographic(const Vector &input, float pole_fade) {
  const Complex coords = stereo(input);
  const float r_sq = coords.re * coords.re + coords.im * coords.im;
  return {coords,
          0,
          0,
          static_cast<uint8_t>(ProjectionBoundary::SINGULAR),
          std::max(0.0f, 1.0f - input.y),
          pole_attenuation(r_sq, pole_fade),
          0};
}

__attribute__((always_inline)) inline ProjectionSample
folded_sinusoidal(const Vector &input, float central_meridian,
                  float pole_fade) {
  const Complex coords =
      projections::folded_sinusoidal(input, central_meridian);
  const float r_sq = coords.re * coords.re + coords.im * coords.im;
  return {coords, static_cast<uint8_t>(input.z < 0.0f), 0,          0,
          1.0f,   pole_attenuation(r_sq, pole_fade),    FOLDED_FLAG};
}

__attribute__((always_inline)) inline ProjectionSample
gnomonic(const Vector &input, float pole_fade, GnomonicHemisphere hemisphere) {
  float y = input.y;
  if (fabsf(y) < GNOMONIC_AXIS_EPS)
    y = y < 0.0f ? -GNOMONIC_AXIS_EPS : GNOMONIC_AXIS_EPS;
  const Complex coords(input.x / y, input.z / y);
  const float r_sq = coords.re * coords.re + coords.im * coords.im;
  const bool in_domain =
      hemisphere == GnomonicHemisphere::FOLDED ||
      (hemisphere == GnomonicHemisphere::FRONT ? input.y >= 0.0f
                                               : input.y < 0.0f);
  return {
      coords,
      static_cast<uint8_t>(input.y < 0.0f),
      static_cast<uint8_t>(input.y < 0.0f),
      static_cast<uint8_t>(static_cast<uint8_t>(ProjectionBoundary::CUT) |
                           static_cast<uint8_t>(ProjectionBoundary::SINGULAR)),
      fabsf(input.y),
      pole_attenuation(r_sq, pole_fade),
      0,
      0,
      0,
      in_domain ? 1.0f : 0.0f};
}

__attribute__((always_inline)) inline ProjectionSample
from_kernel(const projections::ProjectionKernelResult &result,
            float coordinate_scale) {
  return {{result.coords.re * coordinate_scale,
           result.coords.im * coordinate_scale},
          result.region_id,
          result.component_id,
          result.boundary_flags,
          result.fade_edge_distance * fabsf(coordinate_scale),
          1.0f,
          result.flags,
          result.traits,
          result.edge_class};
}

__attribute__((always_inline)) inline ProjectionSample
bonne(const Vector &input, float central_meridian, float standard_parallel,
      float coordinate_scale) {
  return from_kernel(
      projections::bonne_projection(input, central_meridian, standard_parallel),
      coordinate_scale);
}

__attribute__((always_inline)) inline ProjectionSample
peirce(const Vector &input, float central_meridian, uint8_t layout,
       float layout_scroll, bool edge_distance_required,
       float coordinate_scale) {
  return from_kernel(projections::peirce_projection(input, central_meridian,
                                                    layout, layout_scroll,
                                                    edge_distance_required),
                     coordinate_scale);
}

__attribute__((always_inline)) inline ProjectionSample
peirce_fast_square(const Vector &input, float coordinate_scale) {
  return from_kernel(projections::peirce_projection_fast_square(input),
                     coordinate_scale);
}

template <typename State, typename Binding>
concept FrameProvider = Detail::ProviderFor<State, Binding> &&
                        requires(const typename Binding::FrameState &frame) {
                          {
                            State::conjugate(frame)
                          } -> std::same_as<const Quaternion &>;
                        };

__attribute__((always_inline)) inline ProjectionSample
airocean(const Vector &input, float central_meridian, bool horizontal,
         bool edge_distance_required, float coordinate_scale) {
  return from_kernel(projections::airocean_projection(input, central_meridian,
                                                      horizontal,
                                                      edge_distance_required),
                     coordinate_scale);
}

template <typename State, bool North> struct Bonne : ExactPolicy {
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

  static const Quaternion &frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  static ProjectionSample project(const Vector &input,
                                  const FrameState &frame) {
    const float hemisphere = North ? 1.0f : -1.0f;
    return bonne(input, State::central_meridian(frame),
                 hemisphere * State::standard_parallel(frame),
                 State::coordinate_scale(frame));
  }
};

template <typename State> struct Stereographic : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::pole_fade(frame) } -> std::same_as<float>;
      };

  static const Quaternion &frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  static ProjectionSample project(const Vector &input,
                                  const FrameState &frame) {
    return stereographic(input, State::pole_fade(frame));
  }
};

template <typename State> struct FoldedSinusoidal : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::pole_fade(frame) } -> std::same_as<float>;
      };

  static const Quaternion &frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  static ProjectionSample project(const Vector &input,
                                  const FrameState &frame) {
    return folded_sinusoidal(input, State::central_meridian(frame),
                             State::pole_fade(frame));
  }
};

template <typename State, GnomonicHemisphere Hemisphere>
struct Gnomonic : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::pole_fade(frame) } -> std::same_as<float>;
      };

  static const Quaternion &frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  static ProjectionSample project(const Vector &input,
                                  const FrameState &frame) {
    return gnomonic(input, State::pole_fade(frame), Hemisphere);
  }
};

template <typename State, uint8_t Layout, bool EdgeDistanceRequired>
struct Peirce : ExactPolicy {
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
      };

  static const Quaternion &frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  static ProjectionSample project(const Vector &input,
                                  const FrameState &frame) {
    return peirce(input, State::central_meridian(frame), Layout,
                  State::layout_scroll(frame), EdgeDistanceRequired,
                  State::coordinate_scale(frame));
  }
};

template <typename State> struct PeirceFastSquare : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  static constexpr bool EDGE_DISTANCE_UNCONDITIONAL = true;
  static constexpr bool APPROXIMATE = true;
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::PEIRCE_FAST_SQUARE;
  static constexpr std::array<ApproximationMetric, 3> METRICS{{
      {ApproximationDomain::PROJECTED_COORDINATE,
       ApproximationAggregation::MAXIMUM, 1.2e-3f, "plane units"},
      {ApproximationDomain::PROJECTED_EDGE_DISTANCE,
       ApproximationAggregation::MAXIMUM, 2e-4f, "plane units"},
      {ApproximationDomain::FRAMEBUFFER, ApproximationAggregation::MAXIMUM,
       1.0f, "channel code"},
  }};

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
      };

  static bool preconditions(const FrameState &frame) {
    return State::central_meridian(frame) == 0.0f;
  }

  static ProjectionSample project(const Vector &input,
                                  const FrameState &frame) {
    HS_CHECK(preconditions(frame),
             "Peirce fast-square requires a zero central meridian");
    return peirce_fast_square(input, State::coordinate_scale(frame));
  }

  static const Quaternion &frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }
};

template <typename State> struct PeirceSquare : PeirceFastSquare<State> {
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      PeirceFastSquare<State>::template PROVIDER_VALID<CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::layout_scroll(frame) } -> std::same_as<float>;
      };

  static ProjectionSample project(const Vector &input,
                                  const FrameState &frame) {
    if (PeirceFastSquare<State>::preconditions(frame))
      return peirce_fast_square(input, State::coordinate_scale(frame));
    return peirce(input, State::central_meridian(frame), 1,
                  State::layout_scroll(frame), true,
                  State::coordinate_scale(frame));
  }
};

template <typename State, bool Horizontal, bool EdgeDistanceRequired>
struct Airocean : ExactPolicy {
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

  static const Quaternion &frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  static ProjectionSample project(const Vector &input,
                                  const FrameState &frame) {
    return airocean(input, State::central_meridian(frame), Horizontal,
                    EdgeDistanceRequired, State::coordinate_scale(frame));
  }
};

} // namespace Projection

namespace Warp {

static constexpr uint8_t MAX_POLAR_HARMONIC = 16;

struct FlatEnvelope {};
struct ProjectionWeightEnvelope {};
struct EdgeFadeEnvelope {};
struct Euler1 {};
struct Midpoint2 {};
struct Midpoint4 {};
struct LinearPolar {};
struct LogarithmicPolar {};

template <typename State, typename Binding>
concept PreparedProvider =
    Detail::ProviderFor<State, Binding> &&
    requires(const typename Binding::FrameState &frame) {
      State::prepared(frame);
    } &&
    std::is_lvalue_reference_v<decltype(State::prepared(
        std::declval<const typename Binding::FrameState &>()))>;

template <typename State, typename Binding>
concept ParamsPreparedProvider =
    PreparedProvider<State, Binding> &&
    requires(const typename Binding::FrameState &frame) {
      State::params(frame);
    } &&
    std::is_lvalue_reference_v<decltype(State::params(
        std::declval<const typename Binding::FrameState &>()))>;

__attribute__((always_inline)) inline WarpStepResult
finish_closed_form(const Complex &input, const Complex &output) {
  const Complex delta(output.re - input.re, output.im - input.im);
  const float deformation = sqrtf(delta.re * delta.re + delta.im * delta.im);
  return {output, delta, deformation, deformation};
}

template <typename Prepared>
__attribute__((always_inline)) inline WarpStepResult
affine_frame(const Complex &input, const Prepared &prepared) {
  const float c = prepared.rotation_cos;
  const float s = prepared.rotation_sin;
  const auto &affine = prepared.transform.affine;
  const float rx = c * input.re + s * input.im;
  const float ry = -s * input.re + c * input.im;
  return finish_closed_form(
      input, {rx / affine.scale_x - affine.shear * ry / affine.scale_y -
                  affine.translation_x,
              ry / affine.scale_y - affine.translation_y});
}

template <typename Params, typename Prepared>
__attribute__((always_inline)) inline WarpStepResult
wave_shear(const Complex &input, const Params &params, float phase,
           float amplitude, const Prepared &prepared) {
  if (params.strength == 0.0f)
    return {input, Complex(), 0.0f, 0.0f};
  const float c = prepared.rotation_cos;
  const float s = prepared.rotation_sin;
  const float angle =
      params.frequency * (c * input.re + s * input.im) + TWO_PI_F * phase;
  const float offset = amplitude * sinf(angle);
  const Complex delta(-s * offset, c * offset);
  return {{input.re + delta.re, input.im + delta.im},
          delta,
          fabsf(offset),
          fabsf(offset)};
}

template <typename Params, typename Prepared>
__attribute__((always_inline)) inline Complex
mirror_tile_coords(const Complex &input, const Params &params,
                   const Prepared &prepared) {
  const float c = prepared.rotation_cos;
  const float s = prepared.rotation_sin;
  const float offset_x = prepared.transform.mirror.offset_x;
  const float offset_y = prepared.transform.mirror.offset_y;
  const float x = c * input.re + s * input.im + offset_x;
  const float y = -s * input.re + c * input.im + offset_y;
  const float folded_x =
      params.cell_x * (1.0f - 2.0f * fabsf(wrap_t(x / params.cell_x) - 0.5f));
  const float folded_y =
      params.cell_y * (1.0f - 2.0f * fabsf(wrap_t(y / params.cell_y) - 0.5f));
  return {c * folded_x - s * folded_y, s * folded_x + c * folded_y};
}

template <typename Params, typename Prepared>
__attribute__((always_inline)) inline WarpStepResult
mirror_tile(const Complex &input, const Params &params,
            const Prepared &prepared) {
  return finish_closed_form(input, mirror_tile_coords(input, params, prepared));
}

template <typename Params>
__attribute__((always_inline)) inline WarpStepResult
polar_chart(const Complex &input, const Params &params, float phase,
            bool logarithmic, uint8_t harmonic) {
  const float radius = sqrtf(input.re * input.re + input.im * input.im);
  const float radial =
      logarithmic ? logf(std::max(radius, 1.0f / 4096.0f)) : radius;
  return finish_closed_form(
      input, {params.radial_scale * radial + params.radial_phase,
              static_cast<float>(harmonic) * fast_atan2(input.im, input.re) +
                  params.angular_phase + TWO_PI_F * phase});
}

template <typename Params, typename Prepared>
HS_FLASH_MEMBER inline WarpStepResult
vector_noise(const Complex &input, const Params &params, float amplitude,
             const FastNoiseLite &noise, ::NoiseBasis basis,
             const Prepared &prepared) {
  if (params.strength == 0.0f)
    return {input, Complex(), 0.0f, 0.0f};
  const auto &loop = prepared.transform.noise_loop;
  const Vector q(params.scale * input.re + loop.diagonal,
                 params.scale * input.im + loop.diagonal, loop.z);
  float nx;
  float ny;
  if (basis == ::NoiseBasis::SIMPLEX) {
    const Vector field = sample_simplex_vector(noise, q);
    nx = field.x;
    ny = field.y;
  } else {
    nx = sample_noise_vector_channel(noise, basis, q, 0);
    ny = sample_noise_vector_channel(noise, basis, q, 1);
  }
  const float c = prepared.rotation_cos;
  const float s = prepared.rotation_sin;
  const Complex delta(amplitude * (c * nx - s * ny),
                      amplitude * (s * nx + c * ny));
  const float deformation = sqrtf(delta.re * delta.re + delta.im * delta.im);
  return {{input.re + delta.re, input.im + delta.im},
          delta,
          deformation,
          deformation};
}

struct Identity : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static WarpStepResult
  apply(const Complex &input, const ProjectionSample &, const FrameState &) {
    return {input, Complex(), 0.0f, 0.0f};
  }
};

template <typename State> struct AffineFrame : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      PreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::prepared(frame).rotation_cos } -> std::convertible_to<float>;
        { State::prepared(frame).rotation_sin } -> std::convertible_to<float>;
        State::prepared(frame).transform.affine;
      };

  static WarpStepResult apply(const Complex &input, const ProjectionSample &,
                              const FrameState &frame) {
    return affine_frame(input, State::prepared(frame));
  }
};

template <typename State> struct WaveShear : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ParamsPreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).strength } -> std::convertible_to<float>;
        { State::params(frame).frequency } -> std::convertible_to<float>;
        { State::phase(frame) } -> std::same_as<float>;
        { State::prepared(frame).rotation_cos } -> std::convertible_to<float>;
        { State::prepared(frame).rotation_sin } -> std::convertible_to<float>;
      };

  static WarpStepResult apply(const Complex &input, const ProjectionSample &,
                              const FrameState &frame) {
    const auto &params = State::params(frame);
    return wave_shear(input, params, State::phase(frame), params.strength,
                      State::prepared(frame));
  }
};

template <typename State> struct MirrorTile : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ParamsPreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).cell_x } -> std::convertible_to<float>;
        { State::params(frame).cell_y } -> std::convertible_to<float>;
        { State::prepared(frame).rotation_cos } -> std::convertible_to<float>;
        { State::prepared(frame).rotation_sin } -> std::convertible_to<float>;
        State::prepared(frame).transform.mirror;
      };

  static WarpStepResult apply(const Complex &input, const ProjectionSample &,
                              const FrameState &frame) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();
    const WarpStepResult result =
        mirror_tile(input, State::params(frame), State::prepared(frame));
    Instrumentation::template span<ProfileEvent::MIRROR_TILE>(start);
    return result;
  }
};

template <typename State, typename PolarMode, uint8_t Harmonic>
struct PolarChart : ExactPolicy {
  static_assert(Harmonic >= 1 && Harmonic <= MAX_POLAR_HARMONIC);
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ParamsPreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).radial_scale } -> std::convertible_to<float>;
        { State::params(frame).radial_phase } -> std::convertible_to<float>;
        { State::params(frame).angular_phase } -> std::convertible_to<float>;
        { State::phase(frame) } -> std::same_as<float>;
      };

  static WarpStepResult apply(const Complex &input, const ProjectionSample &,
                              const FrameState &frame) {
    return polar_chart(input, State::params(frame), State::phase(frame),
                       std::is_same_v<PolarMode, LogarithmicPolar>, Harmonic);
  }
};

template <typename State, ::NoiseBasis BasisV, typename Envelope>
struct VectorNoise : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ParamsPreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).strength } -> std::convertible_to<float>;
        { State::params(frame).scale } -> std::convertible_to<float>;
        { State::params(frame).edge_width } -> std::convertible_to<float>;
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::prepared(frame).rotation_cos } -> std::convertible_to<float>;
        { State::prepared(frame).rotation_sin } -> std::convertible_to<float>;
        State::prepared(frame).transform.noise_loop;
      };

  static WarpStepResult apply(const Complex &input,
                              const ProjectionSample &projected,
                              const FrameState &frame) {
    const auto &params = State::params(frame);
    float envelope = 1.0f;
    if constexpr (std::is_same_v<Envelope, ProjectionWeightEnvelope>)
      envelope = projected.value_weight;
    else if constexpr (std::is_same_v<Envelope, EdgeFadeEnvelope>)
      envelope = cubic_kernel(projected.fade_edge_distance / params.edge_width);
    return vector_noise(input, params, params.strength * envelope,
                        State::noise(frame), BasisV, State::prepared(frame));
  }
};

} // namespace Warp

namespace Source {

template <typename State, typename Binding>
concept ParamsProvider =
    Detail::ProviderFor<State, Binding> &&
    requires(const typename Binding::FrameState &frame) {
      State::params(frame);
    } &&
    std::is_lvalue_reference_v<decltype(State::params(
        std::declval<const typename Binding::FrameState &>()))>;

template <typename State, typename Binding>
concept StateProvider =
    ParamsProvider<State, Binding> &&
    requires(const typename Binding::FrameState &frame) {
      State::prepared(frame);
    } &&
    std::is_lvalue_reference_v<decltype(State::prepared(
        std::declval<const typename Binding::FrameState &>()))>;

template <typename Prepared>
HS_FLASH_MEMBER inline float twin_wave(const Complex &input,
                                       const Prepared &prepared) {
  const float rotated =
      input.re * prepared.angle_cos + input.im * prepared.angle_sin;
  return 0.5f * (fast_sinf(input.re + prepared.primary) +
                 fast_sinf(rotated + prepared.primary));
}

template <typename Prepared>
HS_FLASH_MEMBER inline float rings(const Complex &input,
                                   const Prepared &prepared) {
  return fast_sinf(sqrtf(input.re * input.re + input.im * input.im) -
                   prepared.primary);
}

template <typename Prepared>
HS_FLASH_MEMBER inline float spiral(const Complex &input,
                                    const Prepared &prepared) {
  const float radius = sqrtf(input.re * input.re + input.im * input.im);
  const float azimuth = fast_atan2(input.im, input.re);
  return fast_sinf(radius - 3.0f * (azimuth + prepared.angle) -
                   prepared.primary);
}

template <typename Params, typename Prepared>
HS_O3_FN inline float grid(const Complex &input, const Params &params,
                           const Prepared &prepared) {
  const float x = input.re * prepared.angle_cos + input.im * prepared.angle_sin;
  const float y =
      -input.re * prepared.angle_sin + input.im * prepared.angle_cos;
  if (params.pattern_mix == 1.0f)
    return fast_sinf(x + prepared.primary) * fast_cosf(y - prepared.secondary);
  float re = x;
  float im = y;
  if (params.complexity != 0.0f) {
    re += params.complexity * fast_sinf(y + prepared.primary);
    im += params.complexity * fast_cosf(x - prepared.secondary);
  }
  const float coupled = fast_sinf(re) * fast_cosf(im);
  if (params.pattern_mix == 0.0f)
    return coupled;
  const float direct =
      fast_sinf(x + prepared.primary) * fast_cosf(y - prepared.secondary);
  return hs::lerp(coupled, direct, params.pattern_mix);
}

template <typename Params>
HS_FLASH_MEMBER inline float primitive_lattice(const Complex &input,
                                               const Params &params) {
  const float x = wrap_t(params.lattice_cell_scale * input.re + 0.5f) - 0.5f;
  const float y = wrap_t(params.lattice_cell_scale * input.im + 0.5f) - 0.5f;
  const float circle = sqrtf(x * x + y * y) - params.lattice_radius;
  const float bx = fabsf(x) - params.lattice_radius;
  const float by = fabsf(y) - params.lattice_radius;
  const float square = sqrtf(std::max(bx, 0.0f) * std::max(bx, 0.0f) +
                             std::max(by, 0.0f) * std::max(by, 0.0f)) +
                       std::min(std::max(bx, by), 0.0f);
  const float distance = hs::lerp(circle, square, params.lattice_shape_blend);
  return 1.0f - 2.0f * ::smooth_ramp(-params.lattice_softness,
                                     params.lattice_softness, distance);
}

template <typename State> struct TwinWave : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::prepared(frame).angle_cos } -> std::convertible_to<float>;
        { State::prepared(frame).angle_sin } -> std::convertible_to<float>;
        { State::prepared(frame).primary } -> std::convertible_to<float>;
      };

  static float sample(const SourceInput &input, const FrameState &frame) {
    const auto &params = State::params(frame);
    return twin_wave(
        stereo_pattern_args(input.warped.coords, params.pattern_freq),
        State::prepared(frame));
  }
};

template <typename State> struct Rings : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::prepared(frame).primary } -> std::convertible_to<float>;
      };

  static float sample(const SourceInput &input, const FrameState &frame) {
    const auto &params = State::params(frame);
    return rings(stereo_pattern_args(input.warped.coords, params.pattern_freq),
                 State::prepared(frame));
  }
};

template <typename State> struct Spiral : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::prepared(frame).angle } -> std::convertible_to<float>;
        { State::prepared(frame).primary } -> std::convertible_to<float>;
      };

  static float sample(const SourceInput &input, const FrameState &frame) {
    const auto &params = State::params(frame);
    return spiral(stereo_pattern_args(input.warped.coords, params.pattern_freq),
                  State::prepared(frame));
  }
};

template <typename State> struct Grid : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::params(frame).pattern_mix } -> std::convertible_to<float>;
        { State::params(frame).complexity } -> std::convertible_to<float>;
        { State::prepared(frame).angle_cos } -> std::convertible_to<float>;
        { State::prepared(frame).angle_sin } -> std::convertible_to<float>;
        { State::prepared(frame).primary } -> std::convertible_to<float>;
        { State::prepared(frame).secondary } -> std::convertible_to<float>;
      };

  static float sample(const SourceInput &input, const FrameState &frame) {
    const auto &params = State::params(frame);
    return grid(stereo_pattern_args(input.warped.coords, params.pattern_freq),
                params, State::prepared(frame));
  }
};

template <typename State> struct PrimitiveLattice : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ParamsProvider<State, CandidateBinding> &&
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

  static float sample(const SourceInput &input, const FrameState &frame) {
    return primitive_lattice(input.warped.coords, State::params(frame));
  }
};

} // namespace Source

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
    return 1.0f - fabsf(2.0f * value - 1.0f);
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

namespace Color {

inline constexpr float UNIT_OPEN_MAX = 0x1.fffffep-1f;

enum class PaletteMapping : uint8_t {
  CUP = 0,
  BELL = 1,
  LINEAR = 2,
  REVERSE = 3
};

enum class BrightnessEnvelope : uint8_t {
  NONE = 0,
  CUP = 1,
  BELL = 2,
  ASCENDING = 3,
  DESCENDING = 4
};

enum class HueMode : uint8_t { NONE = 0, NOISE = 1, PATH_LENGTH = 2 };

struct HueRotationLutView {
  static constexpr int VALUE_STEPS = 64;
  static constexpr int HUE_STEPS = 16;
  static constexpr size_t SIZE = VALUE_STEPS * HUE_STEPS;

  const Pixel *data;
  bool active;
};

struct HueNoiseLutView {
  static constexpr int FACE_COUNT = 6;
  static constexpr int FACE_STEPS = 24;
  static constexpr int FACE_SIZE = FACE_STEPS * FACE_STEPS;
  static constexpr size_t SIZE = FACE_COUNT * FACE_SIZE;

  const int8_t *data;
  bool active;
};

__attribute__((always_inline)) inline float
palette_mapping_coordinate(float value, PaletteMapping mapping, float frequency,
                           float offset) {
  if (mapping == PaletteMapping::LINEAR && frequency == 1.0f && offset == 0.0f)
    return value;
  const float phase =
      wrap_t(std::min(value, UNIT_OPEN_MAX) * frequency + offset);
  switch (mapping) {
  case PaletteMapping::CUP:
    return fabsf(2.0f * phase - 1.0f);
  case PaletteMapping::BELL:
    return 1.0f - fabsf(2.0f * phase - 1.0f);
  case PaletteMapping::LINEAR:
    return phase;
  case PaletteMapping::REVERSE:
    return 1.0f - phase;
  }
  return phase;
}

__attribute__((always_inline)) inline float
brightness_envelope_gain(float value, BrightnessEnvelope envelope,
                         float depth) {
  float shape = 1.0f;
  switch (envelope) {
  case BrightnessEnvelope::NONE:
    return 1.0f;
  case BrightnessEnvelope::CUP:
    shape = fabsf(2.0f * value - 1.0f);
    break;
  case BrightnessEnvelope::BELL:
    shape = 1.0f - fabsf(2.0f * value - 1.0f);
    break;
  case BrightnessEnvelope::ASCENDING:
    shape = value;
    break;
  case BrightnessEnvelope::DESCENDING:
    shape = 1.0f - value;
    break;
  }
  return 1.0f - depth * (1.0f - shape);
}

__attribute__((always_inline)) inline Pixel
sample_hue_rotation_lut(const HueRotationLutView &view, float value,
                        float amount) {
  const float value_position = value * (HueRotationLutView::VALUE_STEPS - 1);
  const int value_low = static_cast<int>(value_position);
  const int value_high =
      std::min(value_low + 1, HueRotationLutView::VALUE_STEPS - 1);
  const uint16_t value_weight =
      frac_to_q16(value_position - static_cast<float>(value_low));

  const float hue_position = amount * HueRotationLutView::HUE_STEPS;
  int hue_low = static_cast<int>(hue_position);
  if (hue_position < static_cast<float>(hue_low))
    --hue_low;
  const uint16_t hue_weight =
      frac_to_q16(hue_position - static_cast<float>(hue_low));
  const int hue_index_low = hue_low & (HueRotationLutView::HUE_STEPS - 1);
  const int hue_index_high =
      (hue_index_low + 1) & (HueRotationLutView::HUE_STEPS - 1);
  const auto sample_row = [&](int row) {
    const int offset = row * HueRotationLutView::HUE_STEPS;
    return view.data[offset + hue_index_low].lerp16(
        view.data[offset + hue_index_high], hue_weight);
  };
  return sample_row(value_low).lerp16(sample_row(value_high), value_weight);
}

HS_FLASH_INLINE inline float sample_hue_noise_lut(const HueNoiseLutView &view,
                                                  const Vector &direction) {
  const float ax = fabsf(direction.x);
  const float ay = fabsf(direction.y);
  const float az = fabsf(direction.z);
  int face;
  float u;
  float v;
  if (ax >= ay && ax >= az) {
    const float inverse = 1.0f / ax;
    face = direction.x >= 0.0f ? 0 : 1;
    u = (direction.x >= 0.0f ? direction.z : -direction.z) * inverse;
    v = direction.y * inverse;
  } else if (ay >= az) {
    const float inverse = 1.0f / ay;
    face = direction.y >= 0.0f ? 2 : 3;
    u = direction.x * inverse;
    v = (direction.y >= 0.0f ? direction.z : -direction.z) * inverse;
  } else {
    const float inverse = 1.0f / az;
    face = direction.z >= 0.0f ? 4 : 5;
    u = (direction.z >= 0.0f ? direction.x : -direction.x) * inverse;
    v = direction.y * inverse;
  }

  constexpr float SCALE =
      0.5f * static_cast<float>(HueNoiseLutView::FACE_STEPS - 1);
  const float x_position = (u + 1.0f) * SCALE;
  const float y_position = (v + 1.0f) * SCALE;
  const int x_low =
      std::min(static_cast<int>(x_position), HueNoiseLutView::FACE_STEPS - 2);
  const int y_low =
      std::min(static_cast<int>(y_position), HueNoiseLutView::FACE_STEPS - 2);
  const float x_fraction = x_position - x_low;
  const float y_fraction = y_position - y_low;
  const int offset = face * HueNoiseLutView::FACE_SIZE +
                     y_low * HueNoiseLutView::FACE_STEPS + x_low;
  const float row_low =
      hs::lerp(static_cast<float>(view.data[offset]),
               static_cast<float>(view.data[offset + 1]), x_fraction);
  const float row_high = hs::lerp(
      static_cast<float>(view.data[offset + HueNoiseLutView::FACE_STEPS]),
      static_cast<float>(view.data[offset + HueNoiseLutView::FACE_STEPS + 1]),
      x_fraction);
  return hs::lerp(row_low, row_high, y_fraction) * (1.0f / 127.0f);
}

template <typename State> struct GeneratedPalette : ExactPolicy {
  static constexpr bool APPROXIMATE = true;
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS;
  static constexpr std::array<ApproximationMetric, 3> METRICS{{
      {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MAXIMUM,
       5400.0f, "channel code"},
      {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MEAN,
       256.0f, "channel code"},
      {ApproximationDomain::FRAMEBUFFER, ApproximationAggregation::MAXIMUM,
       5400.0f, "channel code"},
  }};

  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame, float value) {
        { State::mapping(frame) } -> std::same_as<PaletteMapping>;
        { State::mapping_frequency(frame) } -> std::same_as<float>;
        { State::mapping_phase(frame) } -> std::same_as<float>;
        { State::oscillation_depth(frame) } -> std::same_as<float>;
        { State::oscillation_phase(frame) } -> std::same_as<float>;
        { State::palette(frame, value) } -> std::same_as<Color4>;
        { State::hue_mode(frame) } -> std::same_as<HueMode>;
        { State::hue_shift_amount(frame) } -> std::same_as<float>;
        { State::hue_rotation(frame) } -> std::same_as<HueRotationLutView>;
        { State::hue_noise(frame) } -> std::same_as<HueNoiseLutView>;
        {
          State::brightness_envelope(frame)
        } -> std::same_as<BrightnessEnvelope>;
        { State::brightness_depth(frame) } -> std::same_as<float>;
        { State::opacity_low(frame) } -> std::same_as<float>;
        { State::opacity_high(frame) } -> std::same_as<float>;
      };

  template <typename FrameState>
  HS_FLASH_INLINE static Color4 apply(const MaterialSample &sample,
                                      const FrameState &frame) {
    const float oscillation =
        State::oscillation_depth(frame) *
        fast_sinf(TWO_PI_F * State::oscillation_phase(frame));
    const float palette_value = palette_mapping_coordinate(
        sample.value, State::mapping(frame), State::mapping_frequency(frame),
        State::mapping_phase(frame) + oscillation);
    const HueRotationLutView rotation = State::hue_rotation(frame);
    Color4 color;
    if (rotation.active && State::hue_mode(frame) == HueMode::NOISE) {
      const float amount =
          State::hue_shift_amount(frame) *
          sample_hue_noise_lut(State::hue_noise(frame), sample.sphere);
      color = Color4(sample_hue_rotation_lut(rotation, palette_value, amount),
                     1.0f);
    } else {
      color = State::palette(frame, palette_value);
      if (rotation.active) {
        const float amount =
            wrap_t(State::hue_shift_amount(frame) * sample.path_length);
        if (amount != 0.0f)
          color.color =
              sample_hue_rotation_lut(rotation, palette_value, amount);
      }
    }
    color.color = color.color *
                  brightness_envelope_gain(sample.value,
                                           State::brightness_envelope(frame),
                                           State::brightness_depth(frame));
    color.alpha *=
        sample.coverage * hs::lerp(State::opacity_low(frame),
                                   State::opacity_high(frame), sample.value);
    return color;
  }
};

} // namespace Color

namespace Stage {

template <typename BindingT, typename OrientationState>
struct OuterCamera : Detail::StageContract<BindingT, StageKind::OUTER_CAMERA,
                                           Vector, Vector, false, ExactPolicy> {
  using Base = Detail::StageContract<BindingT, StageKind::OUTER_CAMERA, Vector,
                                     Vector, false, ExactPolicy>;
  using typename Base::FrameState;

  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<OrientationState, BindingT> &&
      requires(const FrameState &frame) {
        {
          OrientationState::conjugate(frame)
        } -> std::same_as<const Quaternion &>;
      };
  static_assert(PROVIDER_VALID,
                "pullback outer camera: malformed orientation provider");

  __attribute__((always_inline)) static Vector run(const Vector &input,
                                                   const FrameState &frame) {
    return rotate(input, OrientationState::conjugate(frame));
  }
};

template <typename BindingT, typename PreLensSurfaceT, typename LensPolicyT,
          typename PostLensSurfaceT, typename ProjectionPolicyT>
struct SurfaceProject
    : Detail::StageContract<BindingT, StageKind::SURFACE_PROJECT, Vector,
                            ProjectionSample, false, PreLensSurfaceT,
                            LensPolicyT, PostLensSurfaceT, ProjectionPolicyT> {
  using Base =
      Detail::StageContract<BindingT, StageKind::SURFACE_PROJECT, Vector,
                            ProjectionSample, false, PreLensSurfaceT,
                            LensPolicyT, PostLensSurfaceT, ProjectionPolicyT>;
  using typename Base::FrameState;
  using PreLensSurface = PreLensSurfaceT;
  using LensPolicy = LensPolicyT;
  using PostLensSurface = PostLensSurfaceT;
  using ProjectionPolicy = ProjectionPolicyT;
  using Instrumentation = typename BindingT::Instrumentation;

  static constexpr bool EDGE_DISTANCE_UNCONDITIONAL = [] {
    if constexpr (requires { ProjectionPolicy::EDGE_DISTANCE_UNCONDITIONAL; })
      return ProjectionPolicy::EDGE_DISTANCE_UNCONDITIONAL;
    else
      return false;
  }();

  static constexpr bool POLICIES_VALID =
      requires(const Vector &input, const FrameState &frame) {
        { PreLensSurface::apply(input, frame) } -> std::same_as<SurfaceResult>;
        { LensPolicy::apply(input, frame) } -> std::same_as<Vector>;
        { PostLensSurface::apply(input, frame) } -> std::same_as<SurfaceResult>;
        {
          ProjectionPolicy::frame_conjugate(frame)
        } -> std::same_as<const Quaternion &>;
        {
          ProjectionPolicy::project(input, frame)
        } -> std::same_as<ProjectionSample>;
      };
  static_assert(POLICIES_VALID,
                "pullback surface/project: malformed policy callable");

  __attribute__((always_inline)) static ProjectionSample
  run(const Vector &input, const FrameState &frame) {
    const SurfaceResult pre = apply_surface<PreLensSurface>(input, frame);
    const Vector lensed = apply_lens(pre.sphere, frame);
    const SurfaceResult post = apply_surface<PostLensSurface>(lensed, frame);
    const auto projection_start = Instrumentation::mark();
    const Vector local =
        rotate(post.sphere, ProjectionPolicy::frame_conjugate(frame));
    ProjectionSample output = ProjectionPolicy::project(local, frame);
    output.sphere = local;
    output.surface_path_length = pre.path_length + post.path_length;
    Instrumentation::template span<ProfileEvent::PROJECTION>(projection_start);
    return output;
  }

private:
  __attribute__((always_inline)) static Vector
  apply_lens(const Vector &input, const FrameState &frame) {
    if constexpr (std::is_same_v<LensPolicy, Lens::Identity>) {
      return LensPolicy::apply(input, frame);
    } else {
      const auto start = Instrumentation::mark();
      const Vector result = LensPolicy::apply(input, frame);
      Instrumentation::template span<ProfileEvent::LENS>(start);
      return result;
    }
  }

  template <typename SurfacePolicy>
  __attribute__((always_inline)) static SurfaceResult
  apply_surface(const Vector &input, const FrameState &frame) {
    if constexpr (std::is_same_v<SurfacePolicy, Surface::Identity>) {
      return SurfacePolicy::apply(input, frame);
    } else {
      const auto start = Instrumentation::mark();
      const SurfaceResult result = SurfacePolicy::apply(input, frame);
      Instrumentation::template span<ProfileEvent::SURFACE_NOISE>(start);
      return result;
    }
  }
};

template <typename BindingT, typename... WarpPolicies>
struct PlanarWarp
    : Detail::StageContract<BindingT, StageKind::PLANAR_WARP, ProjectionSample,
                            SourceInput, false, WarpPolicies...> {
  using Base =
      Detail::StageContract<BindingT, StageKind::PLANAR_WARP, ProjectionSample,
                            SourceInput, false, WarpPolicies...>;
  using typename Base::FrameState;
  using Policies = std::tuple<WarpPolicies...>;
  using Instrumentation = typename BindingT::Instrumentation;
  template <size_t Index>
  using policy_at = std::tuple_element_t<Index, Policies>;

  static constexpr bool POLICIES_VALID =
      (requires(const Complex &input, const ProjectionSample &projected,
                const FrameState &frame) {
        {
          WarpPolicies::apply(input, projected, frame)
        } -> std::same_as<WarpStepResult>;
      } &&
       ...);
  static_assert(POLICIES_VALID,
                "pullback planar warp: malformed policy callable");

  __attribute__((always_inline)) static SourceInput
  run(const ProjectionSample &projected, const FrameState &frame) {
    const auto start = Instrumentation::mark();
    WarpResult warped{projected.coords, Complex(), 0.0f};
    (apply_one<WarpPolicies>(projected, frame, warped), ...);
    Instrumentation::template span<ProfileEvent::PLANAR_WARP>(start);
    return {projected, warped};
  }

private:
  template <typename Policy>
  __attribute__((always_inline)) static void
  apply_one(const ProjectionSample &projected, const FrameState &frame,
            WarpResult &warped) {
    const WarpStepResult step = Policy::apply(warped.coords, projected, frame);
    warped.coords = step.coords;
    warped.net_delta = warped.net_delta + step.delta;
    warped.path_length += step.path_length;
  }
};

template <typename BindingT, typename SourcePolicyT>
struct Source : Detail::StageContract<BindingT, StageKind::SOURCE, SourceInput,
                                      MaterialInput, false, SourcePolicyT> {
  using Base = Detail::StageContract<BindingT, StageKind::SOURCE, SourceInput,
                                     MaterialInput, false, SourcePolicyT>;
  using typename Base::FrameState;
  using SourcePolicy = SourcePolicyT;
  using Instrumentation = typename BindingT::Instrumentation;

  static constexpr bool POLICY_VALID =
      requires(const SourceInput &input, const FrameState &frame) {
        { SourcePolicy::sample(input, frame) } -> std::same_as<float>;
      };
  static_assert(POLICY_VALID, "pullback source: malformed policy callable");

  __attribute__((always_inline)) static MaterialInput
  run(const SourceInput &input, const FrameState &frame) {
    const auto start = Instrumentation::mark();
    const MaterialInput result{input.projected, input.warped,
                               SourcePolicy::sample(input, frame)};
    Instrumentation::template span<ProfileEvent::SOURCE>(start);
    return result;
  }
};

template <typename BindingT, typename WeightPolicyT, typename TransferPolicyT,
          typename CoveragePolicyT>
struct Material
    : Detail::StageContract<BindingT, StageKind::MATERIAL, MaterialInput,
                            MaterialSample, false, WeightPolicyT,
                            TransferPolicyT, CoveragePolicyT> {
  using Base =
      Detail::StageContract<BindingT, StageKind::MATERIAL, MaterialInput,
                            MaterialSample, false, WeightPolicyT,
                            TransferPolicyT, CoveragePolicyT>;
  using typename Base::FrameState;
  using WeightPolicy = WeightPolicyT;
  using TransferPolicy = TransferPolicyT;
  using CoveragePolicy = CoveragePolicyT;
  using Instrumentation = typename BindingT::Instrumentation;

  static constexpr bool POLICIES_VALID = requires(
      float value, const ProjectionSample &projected, const FrameState &frame) {
    { WeightPolicy::apply(value, projected, frame) } -> std::same_as<float>;
    { TransferPolicy::apply(value, frame) } -> std::same_as<float>;
    { CoveragePolicy::apply(value, projected, frame) } -> std::same_as<float>;
  };
  static_assert(POLICIES_VALID, "pullback material: malformed policy callable");

  __attribute__((always_inline)) static MaterialSample
  run(const MaterialInput &input, const FrameState &frame) {
    const auto start = Instrumentation::mark();
    const float weighted =
        WeightPolicy::apply(input.field, input.projected, frame);
    const float unit = Detail::clamp_unit((weighted + 1.0f) * 0.5f);
    const float value = TransferPolicy::apply(unit, frame);
    const float coverage =
        CoveragePolicy::apply(value, input.projected, frame) *
        input.projected.domain_coverage;
    const MaterialSample result{value, coverage, input.projected.sphere,
                                input.projected.surface_path_length +
                                    input.warped.path_length};
    Instrumentation::template span<ProfileEvent::MATERIAL>(start);
    return result;
  }
};

template <typename BindingT, typename ColorPolicyT>
struct Color : Detail::StageContract<BindingT, StageKind::COLOR, MaterialSample,
                                     Color4, true, ColorPolicyT> {
  using Base = Detail::StageContract<BindingT, StageKind::COLOR, MaterialSample,
                                     Color4, true, ColorPolicyT>;
  using typename Base::FrameState;
  using ColorPolicy = ColorPolicyT;
  using Instrumentation = typename BindingT::Instrumentation;

  static constexpr bool POLICY_VALID =
      requires(const MaterialSample &input, const FrameState &frame) {
        { ColorPolicy::apply(input, frame) } -> std::same_as<Color4>;
      };
  static_assert(POLICY_VALID, "pullback color: malformed policy callable");

  __attribute__((always_inline)) static Color4 run(const MaterialSample &input,
                                                   const FrameState &frame) {
    const auto start = Instrumentation::mark();
    const Color4 result = ColorPolicy::apply(input, frame);
    Instrumentation::template span<ProfileEvent::COLOR>(start);
    return result;
  }
};

template <CodeEmission EmissionV, typename StageImplementationT> struct Placed;

template <typename StageImplementationT>
struct Placed<CodeEmission::INLINE_ONLY, StageImplementationT>
    : StageImplementationT {
  using Input = typename StageImplementationT::Input;
  using Output = typename StageImplementationT::Output;
  using FrameState = typename StageImplementationT::FrameState;
  static constexpr CodeEmission EMISSION = CodeEmission::INLINE_ONLY;

  __attribute__((always_inline)) static Output run(const Input &input,
                                                   const FrameState &frame) {
    return StageImplementationT::run(input, frame);
  }
};

template <typename StageImplementationT>
struct Placed<CodeEmission::OUT_OF_LINE_FLASH, StageImplementationT>
    : StageImplementationT {
  using Input = typename StageImplementationT::Input;
  using Output = typename StageImplementationT::Output;
  using FrameState = typename StageImplementationT::FrameState;
  static constexpr CodeEmission EMISSION = CodeEmission::OUT_OF_LINE_FLASH;

  HS_FLASH_MEMBER static Output run(const Input &input,
                                    const FrameState &frame) {
    return StageImplementationT::run(input, frame);
  }
};

template <typename StageImplementationT>
struct Placed<CodeEmission::OUT_OF_LINE_ITCM, StageImplementationT>
    : StageImplementationT {
  using Input = typename StageImplementationT::Input;
  using Output = typename StageImplementationT::Output;
  using FrameState = typename StageImplementationT::FrameState;
  static constexpr CodeEmission EMISSION = CodeEmission::OUT_OF_LINE_ITCM;

  FASTRUN HS_NOINLINE_NOCLONE static Output run(const Input &input,
                                                const FrameState &frame) {
    return StageImplementationT::run(input, frame);
  }
};

} // namespace Stage

} // namespace Pullback
