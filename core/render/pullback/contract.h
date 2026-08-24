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
#include <new>
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
 * @brief Canonical carriers, the ranked family model, the stage-descriptor
 *        contract, validation and the Pipeline composer.
 */

namespace Pullback {

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

/**
 * @brief Rank-0 carrier: a unit view direction plus the accumulated path.
 * @details Must stay all-float: as a 4-float homogeneous aggregate it is
 * returned in s0-s3 across out-of-line boundaries.
 */
struct SphereSample {
  Vector dir;
  float path_length;
};

/** @brief Exactly what a projection computes: regions, fades, weights. */
struct ProjectionProvenance {
  uint8_t region_id;
  uint8_t component_id;
  uint8_t boundary_flags;
  float fade_edge_distance;
  float value_weight;
  uint8_t flags;
  uint8_t traits = 0;
  uint8_t edge_class = 0;
  float domain_coverage = 1.0f;
};

/** @brief The projection policy protocol: planar coordinates plus provenance. */
struct ProjectionResult {
  Complex coords;
  ProjectionProvenance provenance;
};

/**
 * @brief Rank-1 carrier: the working planar coordinate, the projection's
 *        provenance and the sample point the Project combinator recorded.
 * @details Endomorphisms advance `coords` and `path_length` only; `provenance`
 * and `sphere` are immutable once the crossing writes them.
 */
struct PlaneSample {
  Complex coords;
  ProjectionProvenance provenance;
  Vector sphere;
  float path_length;
};

/**
 * @brief Rank-2 carrier: unit field value plus accumulated coverage.
 * @details `value` in [0, 1] is established by the Sample crossing and is
 * each transfer kernel's own obligation thereafter — Kernel::transfer does not
 * re-clamp; `coverage` is in [0, 1] and non-increasing.
 */
struct FieldSample {
  float value;
  float coverage;
  Vector sphere;
  float path_length;
};

// Rank 3 is Color4, straight alpha: endomorphisms are Color4 -> Color4 and
// final premultiplication stays in Scan::Shader, after the chain.

/** @brief The surface policy protocol: a displaced point plus its step length. */
struct SurfaceResult {
  Vector sphere;
  float path_length;
};

/** @brief The warp policy protocol: stepped coordinates plus the step delta. */
struct WarpStepResult {
  Complex coords;
  Complex delta;
  float path_length;
};

namespace Detail {

template <typename... Ts> struct TypeList {
  static constexpr size_t SIZE = sizeof...(Ts);
};

template <typename A, typename B> struct Concat;
template <typename... As, typename... Bs>
struct Concat<TypeList<As...>, TypeList<Bs...>> {
  using Type = TypeList<As..., Bs...>;
};

template <size_t Index, typename List> struct TypeAt;
template <size_t Index, typename... Ts> struct TypeAt<Index, TypeList<Ts...>> {
  using Type = std::tuple_element_t<Index, std::tuple<Ts...>>;
};

} // namespace Detail

/**
 * @brief The closed, ordered carrier set: the single authority every family
 *        fact derives from — membership, rank, and the runtime slot ABI.
 */
using CarrierList =
    Detail::TypeList<SphereSample, PlaneSample, FieldSample, Color4>;

/** Sentinel rank of a type outside the carrier set; keeps family_of total. */
inline constexpr size_t FOREIGN_FAMILY_RANK = static_cast<size_t>(-1);

namespace Detail {

template <typename T, typename List> struct FamilyRank;
template <typename T, typename... Carriers>
struct FamilyRank<T, TypeList<Carriers...>> {
  static constexpr size_t VALUE = [] {
    constexpr bool MATCHES[] = {std::is_same_v<T, Carriers>...};
    for (size_t index = 0; index < sizeof...(Carriers); ++index)
      if (MATCHES[index])
        return index;
    return FOREIGN_FAMILY_RANK;
  }();
};

} // namespace Detail

/**
 * @brief Family rank of @p T: its CarrierList index — Sphere 0, Plane 1,
 *        Field 2, Color 3. Total: a nonmember yields the sentinel rank.
 */
template <typename T>
inline constexpr size_t family_of = Detail::FamilyRank<T, CarrierList>::VALUE;

/** @brief Membership in the closed carrier set. */
template <typename T>
concept CanonicalCarrier = family_of<T> != FOREIGN_FAMILY_RANK;

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

/** @brief Default approximation metadata: exact, no oracle, no metrics. */
struct ApproximationDefaults {
  static constexpr bool APPROXIMATE = false;
  static constexpr bool NON_FLOATING_FIELDS_EXACT = true;
  static constexpr ApproximationOracleId ORACLE = ApproximationOracleId::NONE;
  static constexpr std::array<ApproximationMetric, 0> METRICS{};
};

/** @brief Prepared type of a stage or policy with no per-frame state. */
struct NoPrepared {};

template <typename Stage> consteval bool has_final_framebuffer_metric() {
  for (const ApproximationMetric &metric : Stage::METRICS)
    if (metric.domain == ApproximationDomain::FRAMEBUFFER)
      return true;
  return false;
}

namespace Detail {

/** @brief Whether @p Policy resolves per-frame prepared state. */
template <typename Policy, typename FrameState>
concept PolicyPrepares = requires(const FrameState &frame) {
  { Policy::prepare(frame) } -> std::same_as<typename Policy::Prepared>;
};

template <typename Policy, typename FrameState> struct PolicyPreparedImpl {
  using Type = NoPrepared;
};

template <typename Policy, typename FrameState>
  requires PolicyPrepares<Policy, FrameState>
struct PolicyPreparedImpl<Policy, FrameState> {
  using Type = typename Policy::Prepared;
};

/** @brief The policy's Prepared type, or NoPrepared when it declares none. */
template <typename Policy, typename FrameState>
using PolicyPrepared = typename PolicyPreparedImpl<Policy, FrameState>::Type;

/** @brief Resolves one policy's prepared state; empty when it declares none. */
template <typename Policy, typename FrameState>
__attribute__((always_inline)) inline PolicyPrepared<Policy, FrameState>
prepare_policy(const FrameState &frame) {
  if constexpr (PolicyPrepares<Policy, FrameState>)
    return Policy::prepare(frame);
  else
    return {};
}

template <typename... Policies> struct FirstApproximate {
  using Type = ApproximationDefaults;
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

template <typename Tuple> struct TupleApproximation;
template <typename... Policies>
struct TupleApproximation<std::tuple<Policies...>>
    : CombinedApproximation<Policies...> {};

template <typename Policy, typename Binding>
consteval bool policy_provider_valid() {
  if constexpr (requires { Policy::template PROVIDER_VALID<Binding>; })
    return Policy::template PROVIDER_VALID<Binding>;
  else
    return true;
}

template <typename Tuple, typename Binding> struct PoliciesBindable;
template <typename... Policies, typename Binding>
struct PoliciesBindable<std::tuple<Policies...>, Binding>
    : std::bool_constant<(policy_provider_valid<Policies, Binding>() && ...)> {
};

template <typename T> struct IsPolicyTuple : std::false_type {};
template <typename... Policies>
struct IsPolicyTuple<std::tuple<Policies...>> : std::true_type {};

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

/** Probe binding for naming a descriptor's Bind without instantiating it. */
struct ProbeBinding {
  struct FrameState {};
  using Instrumentation = NoInstrumentation;
};

} // namespace Detail

/** @brief Shape of the public stage-descriptor contract. */
template <typename Descriptor>
concept StageDescriptor = requires {
  typename Descriptor::Input;
  typename Descriptor::Output;
  typename Descriptor::Policies;
  typename Descriptor::template Bind<Detail::ProbeBinding>;
} && Detail::IsPolicyTuple<typename Descriptor::Policies>::value;

/**
 * @brief Whether @p Descriptor's policies and providers agree with @p Binding.
 * @details Non-asserting: evaluated before the bind step instantiates
 * anything, so the pipeline's named BINDINGS assertion reports a foreign
 * binding instead of a template-formation abort inside Bind.
 */
template <typename Descriptor, typename Binding>
consteval bool descriptor_bindable() {
  bool valid =
      Detail::PoliciesBindable<typename Descriptor::Policies, Binding>::value;
  if constexpr (requires { Descriptor::template PROVIDER_VALID<Binding>; })
    valid = valid && Descriptor::template PROVIDER_VALID<Binding>;
  return valid;
}

namespace Detail {

template <typename Descriptor, typename Binding>
concept DescriptorPrepares =
    requires(const typename Binding::FrameState &frame) {
      Descriptor::template prepare<Binding>(frame);
    };

template <typename Descriptor, typename Binding, bool Has>
struct DescriptorPreparedImpl {
  using Type = NoPrepared;
};

template <typename Descriptor, typename Binding>
struct DescriptorPreparedImpl<Descriptor, Binding, true> {
  using Type = decltype(Descriptor::template prepare<Binding>(
      std::declval<const typename Binding::FrameState &>()));
};

/** @brief The descriptor's bound Prepared type: prepare()'s return, or none. */
template <typename Descriptor, typename Binding>
using DescriptorPrepared = typename DescriptorPreparedImpl<
    Descriptor, Binding, DescriptorPrepares<Descriptor, Binding>>::Type;

} // namespace Detail

namespace Stage {

/**
 * @brief Helper base every stage descriptor derives from.
 * @details Derives the carrier pair and forwards Bind into the descriptor's
 * binding-templated `run` (and optional `prepare`) statics — which is how an
 * unbound descriptor defines execution that needs the binding-dependent
 * FrameState and Instrumentation.
 * @tparam Derived The descriptor deriving from this base.
 * @tparam InputT The stage's input carrier.
 * @tparam OutputT The stage's output carrier.
 */
template <typename Derived, typename InputT, typename OutputT> struct Contract {
  using Input = InputT;
  using Output = OutputT;

  template <typename BindingT>
  struct Bind : Detail::TupleApproximation<typename Derived::Policies> {
    using Descriptor = Derived;
    using Binding = BindingT;
    using FrameState = typename BindingT::FrameState;
    using Input = InputT;
    using Output = OutputT;
    using Prepared = Detail::DescriptorPrepared<Derived, BindingT>;

    // Backstop only: pipeline assembly checks descriptor_bindable() first and
    // reports through the named BINDINGS assertion.
    static_assert(descriptor_bindable<Derived, BindingT>(),
                  "pullback stage: malformed or foreign provider");

    HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
      if constexpr (Detail::DescriptorPrepares<Derived, BindingT>)
        return Derived::template prepare<BindingT>(frame);
      else
        return {};
    }

    __attribute__((always_inline)) static Output
    run(const Input &input, const FrameState &frame, const Prepared &prepared) {
      return Derived::template run<BindingT>(input, frame, prepared);
    }
  };
};

/**
 * @brief Placement node: a contiguous run of stages emitted as one call unit
 *        under the given emission.
 * @details Invisible to the semantic leaf view — adding or removing a
 * placement wrapper never changes whether a chain validates or what a
 * predicate matches; placement affects exactly prepared-state layout and code
 * emission. Placements are deliberate, ITCM-ledger-scored decisions, always
 * written by the author.
 */
template <CodeEmission EmissionV, typename... Stages> struct Placed {
  static constexpr CodeEmission EMISSION = EmissionV;
};

} // namespace Stage

namespace Detail {

template <typename... Ts> struct FilterVoid {
  using Type = TypeList<>;
};

template <typename Head, typename... Tail> struct FilterVoid<Head, Tail...> {
  using Rest = typename FilterVoid<Tail...>::Type;
  using Type = std::conditional_t<std::is_void_v<Head>, Rest,
                                  typename Concat<TypeList<Head>, Rest>::Type>;
};

template <CodeEmission EmissionV, typename List> struct RepackPlaced;
template <CodeEmission EmissionV, typename... Stages>
struct RepackPlaced<EmissionV, TypeList<Stages...>> {
  using Type = Stage::Placed<EmissionV, Stages...>;
};

/** Filters void entries; an empty Placed collapses to void like its members. */
template <typename Entry> struct NormalizeNode {
  using Type = Entry;
};

template <CodeEmission EmissionV, typename... Children>
struct NormalizeNode<Stage::Placed<EmissionV, Children...>> {
  using Filtered =
      typename FilterVoid<typename NormalizeNode<Children>::Type...>::Type;
  using Type =
      std::conditional_t<Filtered::SIZE == 0, void,
                         typename RepackPlaced<EmissionV, Filtered>::Type>;
};

template <typename... Entries> struct NormalizeEntries {
  using Type =
      typename FilterVoid<typename NormalizeNode<Entries>::Type...>::Type;
};

template <typename... Lists> struct ConcatAll {
  using Type = TypeList<>;
};

template <typename Head, typename... Tail> struct ConcatAll<Head, Tail...> {
  using Type = typename Concat<Head, typename ConcatAll<Tail...>::Type>::Type;
};

/** Semantic leaf projection: placement nodes flatten to their children. */
template <typename Node> struct NodeLeaves {
  using Type = TypeList<Node>;
};

template <CodeEmission EmissionV, typename... Children>
struct NodeLeaves<Stage::Placed<EmissionV, Children...>> {
  using Type = typename ConcatAll<typename NodeLeaves<Children>::Type...>::Type;
};

template <typename NodeList> struct Leaves;
template <typename... Nodes> struct Leaves<TypeList<Nodes...>> {
  using Type = typename ConcatAll<typename NodeLeaves<Nodes>::Type...>::Type;
};

/** Bound-stage placeholder when descriptor_bindable() fails; keeps the
    pipeline compilable so the named BINDINGS assertion is what reports. */
template <typename Descriptor, typename Binding> struct UnboundStage {
  using Input = typename Descriptor::Input;
  using Output = typename Descriptor::Output;
  using FrameState = typename Binding::FrameState;
  using Prepared = NoPrepared;

  static constexpr Prepared prepare(const FrameState &) { return {}; }
};

template <typename Node, typename Binding> struct BindNode {
  using Type = std::conditional_t<descriptor_bindable<Node, Binding>(),
                                  typename Node::template Bind<Binding>,
                                  UnboundStage<Node, Binding>>;
};

template <CodeEmission EmissionV, typename Binding, typename... Children>
struct BoundPlacedBase {
  using BoundTuple = std::tuple<typename BindNode<Children, Binding>::Type...>;
  static constexpr size_t CHILD_COUNT = sizeof...(Children);
  template <size_t Index>
  using child_at = std::tuple_element_t<Index, BoundTuple>;
  using Input = typename child_at<0>::Input;
  using Output = typename child_at<CHILD_COUNT - 1>::Output;
  using FrameState = typename Binding::FrameState;
  /** Nested per placement node: an out-of-line call receives one contiguous
      prepared sub-tuple. */
  using Prepared =
      std::tuple<typename BindNode<Children, Binding>::Type::Prepared...>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return Prepared{BindNode<Children, Binding>::Type::prepare(frame)...};
  }

protected:
  template <size_t Index, typename In>
  __attribute__((always_inline)) static Output
  run_children(const In &input, const FrameState &frame,
               const Prepared &prepared) {
    using Current = child_at<Index>;
    static_assert(std::is_same_v<In, typename Current::Input>,
                  "pullback pipeline: stage input mismatch");
    const typename Current::Output output =
        Current::run(input, frame, std::get<Index>(prepared));
    if constexpr (Index + 1 == CHILD_COUNT)
      return output;
    else
      return run_children<Index + 1>(output, frame, prepared);
  }
};

template <CodeEmission EmissionV, typename Binding, typename... Children>
struct BoundPlaced;

template <typename Binding, typename... Children>
struct BoundPlaced<CodeEmission::INLINE_ONLY, Binding, Children...>
    : BoundPlacedBase<CodeEmission::INLINE_ONLY, Binding, Children...> {
  using Base = BoundPlacedBase<CodeEmission::INLINE_ONLY, Binding, Children...>;
  static constexpr CodeEmission EMISSION = CodeEmission::INLINE_ONLY;

  __attribute__((always_inline)) static typename Base::Output
  run(const typename Base::Input &input, const typename Base::FrameState &frame,
      const typename Base::Prepared &prepared) {
    return Base::template run_children<0>(input, frame, prepared);
  }
};

template <typename Binding, typename... Children>
struct BoundPlaced<CodeEmission::OUT_OF_LINE_FLASH, Binding, Children...>
    : BoundPlacedBase<CodeEmission::OUT_OF_LINE_FLASH, Binding, Children...> {
  using Base =
      BoundPlacedBase<CodeEmission::OUT_OF_LINE_FLASH, Binding, Children...>;
  static constexpr CodeEmission EMISSION = CodeEmission::OUT_OF_LINE_FLASH;

  __attribute__((noinline)) HS_FLASH_MEMBER static typename Base::Output
  run(const typename Base::Input &input, const typename Base::FrameState &frame,
      const typename Base::Prepared &prepared) {
    return Base::template run_children<0>(input, frame, prepared);
  }
};

template <typename Binding, typename... Children>
struct BoundPlaced<CodeEmission::OUT_OF_LINE_ITCM, Binding, Children...>
    : BoundPlacedBase<CodeEmission::OUT_OF_LINE_ITCM, Binding, Children...> {
  using Base =
      BoundPlacedBase<CodeEmission::OUT_OF_LINE_ITCM, Binding, Children...>;
  static constexpr CodeEmission EMISSION = CodeEmission::OUT_OF_LINE_ITCM;

  FASTRUN HS_NOINLINE_NOCLONE static typename Base::Output
  run(const typename Base::Input &input, const typename Base::FrameState &frame,
      const typename Base::Prepared &prepared) {
    return Base::template run_children<0>(input, frame, prepared);
  }
};

template <CodeEmission EmissionV, typename... Children, typename Binding>
struct BindNode<Stage::Placed<EmissionV, Children...>, Binding> {
  using Type = BoundPlaced<EmissionV, Binding, Children...>;
};

} // namespace Detail

template <typename Binding, typename... Stages>
consteval bool binding_extra_validation() {
  if constexpr (requires {
                  Binding::template ExtraValidation<Stages...>::value;
                })
    return Binding::template ExtraValidation<Stages...>::value;
  else
    return true;
}

namespace Detail {

template <typename LeafList> struct LeafShape;
template <typename... Ls> struct LeafShape<TypeList<Ls...>> {
  static constexpr size_t COUNT = sizeof...(Ls);
  static constexpr bool NONEMPTY = COUNT > 0;
  static constexpr bool CONTRACTS = (StageDescriptor<Ls> && ...);
};

template <bool ShapeValid, typename LeafList> struct LeafCarrierRows {
  static constexpr bool CANONICAL = false;
  static constexpr bool MONOTONE = false;
  static constexpr bool CARRIERS = false;
  static constexpr bool ENTRY = false;
  static constexpr bool EXIT = false;
};

template <typename... Ls> struct LeafCarrierRows<true, TypeList<Ls...>> {
private:
  using Tuple = std::tuple<Ls...>;

  template <size_t... Is>
  static consteval bool adjacent(std::index_sequence<Is...>) {
    return (
        std::is_same_v<typename std::tuple_element_t<Is, Tuple>::Output,
                       typename std::tuple_element_t<Is + 1, Tuple>::Input> &&
        ...);
  }

public:
  static constexpr bool CANONICAL = ((CanonicalCarrier<typename Ls::Input> &&
                                      CanonicalCarrier<typename Ls::Output>) &&
                                     ...);
  static constexpr bool MONOTONE =
      ((family_of<typename Ls::Input> <= family_of<typename Ls::Output>) &&
       ...);
  static constexpr bool CARRIERS =
      adjacent(std::make_index_sequence<sizeof...(Ls) - 1>{});
  static constexpr bool ENTRY =
      std::is_same_v<typename std::tuple_element_t<0, Tuple>::Input,
                     SphereSample>;
  static constexpr bool EXIT = std::is_same_v<
      typename std::tuple_element_t<sizeof...(Ls) - 1, Tuple>::Output, Color4>;
};

template <bool ShapeValid, typename Binding, typename LeafList>
struct LeafBindingRows {
  static constexpr bool BINDINGS = false;
  static constexpr bool EMPTY_DESCRIPTORS = false;
  static constexpr bool EXTRA_VALIDATION = false;
};

template <typename Binding, typename... Ls>
struct LeafBindingRows<true, Binding, TypeList<Ls...>> {
  static constexpr bool BINDINGS = (descriptor_bindable<Ls, Binding>() && ...);
  static constexpr bool EMPTY_DESCRIPTORS = (std::is_empty_v<Ls> && ...);
  static constexpr bool EXTRA_VALIDATION =
      binding_extra_validation<Binding, Ls...>();
};

template <typename Bound> consteval bool bound_run_returns() {
  return requires(const typename Bound::Input &input,
                  const typename Bound::FrameState &frame,
                  const typename Bound::Prepared &prepared) {
    {
      Bound::run(input, frame, prepared)
    } -> std::same_as<typename Bound::Output>;
  };
}

template <typename Bound> consteval bool bound_prepares() {
  return requires(const typename Bound::FrameState &frame) {
    { Bound::prepare(frame) } -> std::same_as<typename Bound::Prepared>;
  };
}

template <typename Bound> consteval bool bound_approximation_wellformed() {
  return (Bound::APPROXIMATE && Bound::ORACLE != ApproximationOracleId::NONE &&
          Bound::METRICS.size() != 0 && Bound::NON_FLOATING_FIELDS_EXACT &&
          has_final_framebuffer_metric<Bound>()) ||
         (!Bound::APPROXIMATE && Bound::ORACLE == ApproximationOracleId::NONE &&
          Bound::METRICS.size() == 0);
}

template <bool BindingsValid, typename Binding, typename LeafList>
struct LeafCallableRows {
  static constexpr bool RUN_RETURNS = false;
  static constexpr bool PREPARES = false;
  static constexpr bool APPROXIMATIONS = false;
};

template <typename Binding, typename... Ls>
struct LeafCallableRows<true, Binding, TypeList<Ls...>> {
  static constexpr bool RUN_RETURNS =
      (bound_run_returns<typename Ls::template Bind<Binding>>() && ...);
  static constexpr bool PREPARES =
      (bound_prepares<typename Ls::template Bind<Binding>>() && ...);
  static constexpr bool APPROXIMATIONS =
      (bound_approximation_wellformed<typename Ls::template Bind<Binding>>() &&
       ...);
};

} // namespace Detail

/**
 * @brief The staged, named-boolean validation surface over a chain's
 *        flattened semantic leaf list.
 * @details Earlier levels gate later ones: descriptor contract shape, then
 * CANONICAL, then the rank rows, then the bound callable checks — so a
 * malformed descriptor reports through its own named row instead of
 * detonating a later fold. Placement wrappers are invisible to every row.
 */
template <typename Binding, typename... Entries> struct PipelineValidation {
  using NodeList = typename Detail::NormalizeEntries<Entries...>::Type;
  using LeafList = typename Detail::Leaves<NodeList>::Type;

private:
  using Shape = Detail::LeafShape<LeafList>;
  static constexpr bool SHAPE = Shape::NONEMPTY && Shape::CONTRACTS;
  using CarrierRows = Detail::LeafCarrierRows<SHAPE, LeafList>;
  using BindingRows = Detail::LeafBindingRows<SHAPE, Binding, LeafList>;
  using CallableRows = Detail::LeafCallableRows<SHAPE && BindingRows::BINDINGS,
                                                Binding, LeafList>;

public:
  static constexpr size_t LEAF_COUNT = Shape::COUNT;
  static constexpr bool NONEMPTY = Shape::NONEMPTY;
  static constexpr bool CONTRACTS = Shape::CONTRACTS;
  static constexpr bool CANONICAL = CarrierRows::CANONICAL;
  static constexpr bool MONOTONE = CarrierRows::MONOTONE;
  static constexpr bool CARRIERS = CarrierRows::CARRIERS;
  static constexpr bool ENTRY = CarrierRows::ENTRY;
  static constexpr bool EXIT = CarrierRows::EXIT;
  static constexpr bool BINDINGS = BindingRows::BINDINGS;
  static constexpr bool EMPTY_DESCRIPTORS = BindingRows::EMPTY_DESCRIPTORS;
  static constexpr bool EXTRA_VALIDATION = BindingRows::EXTRA_VALIDATION;
  static constexpr bool RUN_RETURNS = CallableRows::RUN_RETURNS;
  static constexpr bool PREPARES = CallableRows::PREPARES;
  static constexpr bool APPROXIMATIONS = CallableRows::APPROXIMATIONS;
};

template <typename Stage, typename Key>
concept StageMatchesKey = requires(const Key &key) {
  { Stage::implements(key) } -> std::convertible_to<bool>;
};

namespace Detail {

template <template <typename> class Predicate, typename... Ls>
struct FirstMatching {
  using Type = void;
};

template <template <typename> class Predicate, typename Head, typename... Tail>
struct FirstMatching<Predicate, Head, Tail...> {
  using Type =
      std::conditional_t<Predicate<Head>::value, Head,
                         typename FirstMatching<Predicate, Tail...>::Type>;
};

template <typename LeafList> struct LeafOps;
template <typename... Ls> struct LeafOps<TypeList<Ls...>> {
  template <typename Key>
    requires(StageMatchesKey<Ls, Key> && ...)
  static constexpr bool implements(const Key &key) {
    return (Ls::implements(key) && ...);
  }

  template <template <typename> class Predicate>
  static constexpr bool ANY = (Predicate<Ls>::value || ...);

  template <template <typename> class Predicate>
  using Matching = typename FirstMatching<Predicate, Ls...>::Type;
};

template <bool Valid, typename Binding, typename NodeList> struct PipelineCore {
  using PreparedTuple = std::tuple<>;
  static constexpr size_t NODE_COUNT = 0;
  template <size_t Index> using NodeAt = void;

  static PreparedTuple prepare_stages(const typename Binding::FrameState &) {
    return {};
  }
};

template <typename Binding, typename... Nodes>
struct PipelineCore<true, Binding, TypeList<Nodes...>> {
  using FrameState = typename Binding::FrameState;
  using BoundNodes = std::tuple<typename BindNode<Nodes, Binding>::Type...>;
  static constexpr size_t NODE_COUNT = sizeof...(Nodes);
  template <size_t Index>
  using NodeAt = std::tuple_element_t<Index, BoundNodes>;
  /** Per-node prepared state; a placement node's entry nests its children. */
  using PreparedTuple =
      std::tuple<typename BindNode<Nodes, Binding>::Type::Prepared...>;

  HS_FLASH_MEMBER static PreparedTuple prepare_stages(const FrameState &ctx) {
    return PreparedTuple{BindNode<Nodes, Binding>::Type::prepare(ctx)...};
  }

  template <size_t Index, typename Input>
  __attribute__((always_inline)) static Color4
  run_stage(const Input &input, const FrameState &ctx,
            const PreparedTuple &prepared) {
    using CurrentNode = NodeAt<Index>;
    static_assert(std::is_same_v<Input, typename CurrentNode::Input>,
                  "pullback pipeline: stage input mismatch");
    const typename CurrentNode::Output output =
        CurrentNode::run(input, ctx, std::get<Index>(prepared));
    if constexpr (Index + 1 == NODE_COUNT)
      return output;
    else
      return run_stage<Index + 1>(output, ctx, prepared);
  }
};

} // namespace Detail

/**
 * @brief The ranked pullback pipeline: one binding, then a chain of stage
 *        descriptors and placement nodes.
 * @details Binding appears exactly once — the first parameter binds the whole
 * list. `void` entries and empty placement groups vanish, which is the entire
 * hook a derivation layer needs; concrete pipelines never contain a
 * conditional. The chain must be non-decreasing in family rank, adjacent
 * carriers must agree, the first stage consumes SphereSample and the last
 * produces Color4. The semantic leaf view (validation, predicates,
 * STAGE_COUNT, stage_at) and the structural placement view (prepare/run over
 * execution nodes) share no indices.
 */
template <typename BindingT, typename... Entries> struct Pipeline {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  using Validation = PipelineValidation<BindingT, Entries...>;
  using NodeList = typename Validation::NodeList;
  using LeafList = typename Validation::LeafList;

  static_assert(Validation::NONEMPTY,
                "pullback pipeline: at least one stage is required");
  static_assert(!Validation::NONEMPTY || Validation::CONTRACTS,
                "pullback pipeline: missing or mistyped stage descriptor");
  static_assert(!Validation::NONEMPTY || !Validation::CONTRACTS ||
                    Validation::CANONICAL,
                "pullback pipeline: stage carrier is outside the canonical "
                "set (SphereSample, PlaneSample, FieldSample, Color4)");
  static_assert(!Validation::NONEMPTY || !Validation::CONTRACTS ||
                    !Validation::CANONICAL || Validation::MONOTONE,
                "pullback pipeline: a stage may not decrease its family "
                "rank — families are ordered Sphere, Plane, Field, Color");
  static_assert(!Validation::NONEMPTY || !Validation::CONTRACTS ||
                    !Validation::CANONICAL || Validation::CARRIERS,
                "pullback pipeline: carrier mismatch");
  static_assert(!Validation::NONEMPTY || !Validation::CONTRACTS ||
                    !Validation::CANONICAL || Validation::ENTRY,
                "pullback pipeline: the first stage must consume SphereSample");
  static_assert(!Validation::NONEMPTY || !Validation::CONTRACTS ||
                    !Validation::CANONICAL || Validation::EXIT,
                "pullback pipeline: the last stage must produce Color4");
  static_assert(!Validation::NONEMPTY || !Validation::CONTRACTS ||
                    Validation::BINDINGS,
                "pullback pipeline: stage binding mismatch");
  static_assert(!Validation::NONEMPTY || !Validation::CONTRACTS ||
                    Validation::EMPTY_DESCRIPTORS,
                "pullback pipeline: stage descriptors must be stateless");
  static_assert(!Validation::NONEMPTY || !Validation::CONTRACTS ||
                    !Validation::BINDINGS || Validation::RUN_RETURNS,
                "pullback pipeline: wrong stage return type");
  static_assert(!Validation::NONEMPTY || !Validation::CONTRACTS ||
                    !Validation::BINDINGS || Validation::PREPARES,
                "pullback pipeline: wrong stage prepare() return type");
  static_assert(!Validation::NONEMPTY || !Validation::CONTRACTS ||
                    !Validation::BINDINGS || Validation::APPROXIMATIONS,
                "pullback pipeline: malformed approximation metadata");
  static_assert(!Validation::NONEMPTY || !Validation::CONTRACTS ||
                    Validation::EXTRA_VALIDATION,
                "pullback pipeline: consumer validation failed");

  /** Semantic leaf view: counts and indexes leaves, never execution nodes. */
  static constexpr size_t STAGE_COUNT = Validation::LEAF_COUNT;
  template <size_t Index>
  using stage_at = typename Detail::TypeAt<Index, LeafList>::Type;

private:
  static constexpr bool CORE_VALID =
      Validation::NONEMPTY && Validation::CONTRACTS;
  using Core = Detail::PipelineCore<CORE_VALID, BindingT, NodeList>;

public:
  /** Execution view: prepare and run recurse over placement nodes. */
  static constexpr size_t NODE_COUNT = Core::NODE_COUNT;
  template <size_t Index> using node_at = typename Core::template NodeAt<Index>;

  /** Per-node prepared state, nested per placement group. */
  using PreparedTuple = typename Core::PreparedTuple;
  static_assert(std::is_trivially_destructible_v<PreparedTuple>,
                "pullback pipeline: prepared state must be trivially "
                "destructible");

  /**
   * @brief One frame's complete shading state: the public frame context plus
   *        the stages' private prepared state.
   */
  struct Frame {
    FrameState ctx;
    PreparedTuple prepared;
  };

  /** @brief Resolves every stage's prepared state from @p ctx. */
  HS_FLASH_MEMBER static PreparedTuple prepare_stages(const FrameState &ctx) {
    return Core::prepare_stages(ctx);
  }

  /** @brief Bundles @p ctx with the stages' prepared state. */
  HS_FLASH_MEMBER static Frame prepare(const FrameState &ctx) {
    return {ctx, prepare_stages(ctx)};
  }

  /** @brief Type-erased prepare for dynamic program dispatch. */
  HS_FLASH_MEMBER static void prepare_into(const FrameState &ctx,
                                           void *storage) {
    new (storage) PreparedTuple{prepare_stages(ctx)};
  }

  /** @brief Seeds the entry carrier and runs the chain over @p view. */
  __attribute__((always_inline)) static Color4
  evaluate(const Vector &view, const FrameState &ctx,
           const PreparedTuple &prepared) {
    return Core::template run_stage<0>(SphereSample{view, 0.0f}, ctx, prepared);
  }

  HS_FLASH_MEMBER static Color4 shade(const Vector &view, const Frame &frame) {
    return evaluate(view, frame.ctx, frame.prepared);
  }

  /** @brief Type-erased shade over prepare_into()'s storage. */
  HS_FLASH_MEMBER static Color4 shade_prepared(const Vector &view,
                                               const FrameState &ctx,
                                               const void *storage) {
    return evaluate(view, ctx, *static_cast<const PreparedTuple *>(storage));
  }

  template <typename Key> static constexpr bool implements(const Key &key) {
    return Detail::LeafOps<LeafList>::implements(key);
  }

  /** @brief Whether any leaf stage satisfies @p Predicate. */
  template <template <typename> class Predicate>
  static constexpr bool any_stage =
      Detail::LeafOps<LeafList>::template ANY<Predicate>;

  /** @brief The first leaf stage satisfying @p Predicate, or void. */
  template <template <typename> class Predicate>
  using stage_matching =
      typename Detail::LeafOps<LeafList>::template Matching<Predicate>;
};

} // namespace Pullback
