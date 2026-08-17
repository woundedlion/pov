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
#include <span>
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
 * @file contract.h
 * @brief Stage contracts, carriers, validation and the Pipeline composer.
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

using ProjectionBoundary = projections::ProjectionBoundary;

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
  float path_length;
};

struct WarpResult {
  Complex coords;
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

__attribute__((always_inline)) inline float clamp_unit(float value) {
  if (value <= 0.0f)
    return 0.0f;
  if (value >= 1.0f)
    return 1.0f;
  return value;
}

__attribute__((always_inline)) inline float smooth_ramp(float low, float high,
                                                        float value) {
  if (low == high)
    return static_cast<float>(value > low);
  const float t = clamp_unit((value - low) / (high - low));
  return t * t * (3.0f - 2.0f * t);
}

} // namespace Detail

} // namespace Pullback
