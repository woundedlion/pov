/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file composed_resources.h
 * @brief Per-key resources for composed pullback stages. */

#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

#include "math/noise_field.h"
#include "render/pullback/source.h"
#include "render/pullback/stage.h"
#include "render/pullback/warp.h"

namespace Pullback {

/** @brief Compile-time identity of one composed parameter/resource instance. */
template <size_t N> struct ResourceKey {
  char text[N];
  constexpr ResourceKey(const char (&value)[N]) {
    for (size_t i = 0; i < N; ++i)
      text[i] = value[i];
  }
  constexpr std::string_view view() const { return {text, N - 1}; }
};

enum class ResourceKind {
  SOURCE,
  PROJECTION,
  WARP,
  SURFACE,
  LENS,
  VALUE,
  COLOR
};

/** @brief The storage required by one provider instance. */
template <ResourceKey Key, typename FamilyT, ResourceKind KindV>
struct ParameterResource {
  static constexpr auto KEY = Key;
  static constexpr ResourceKind KIND = KindV;
  using Family = FamilyT;
  static constexpr size_t ORDER = [] {
    constexpr std::string_view INSTANCE_ORDER[] = {
        "source",  "projection", "outer_warp", "inner_warp",
        "surface", "lens",       "value",      "color"};
    for (size_t i = 0; i < std::size(INSTANCE_ORDER); ++i)
      if (Key.view() == INSTANCE_ORDER[i])
        return i;
    return std::size(INSTANCE_ORDER);
  }();
};

namespace ComposedDetail {

template <typename... T> struct ResourceList {};

template <typename List, typename Resource> struct InsertResource;
template <typename Resource> struct InsertResource<ResourceList<>, Resource> {
  using Type = ResourceList<Resource>;
};
template <typename First, typename... Rest, typename Resource>
struct InsertResource<ResourceList<First, Rest...>, Resource> {
  static constexpr bool SAME_KEY = First::KEY.view() == Resource::KEY.view();
  static_assert(!SAME_KEY || std::is_same_v<First, Resource>,
                "composed resource key has conflicting families or roles");
  template <typename... T>
  static ResourceList<First, T...> prepend(ResourceList<T...>);
  using Type = std::conditional_t<
      SAME_KEY, ResourceList<First, Rest...>,
      std::conditional_t<(Resource::ORDER < First::ORDER),
                         ResourceList<Resource, First, Rest...>,
                         decltype(prepend(
                             typename InsertResource<ResourceList<Rest...>,
                                                     Resource>::Type{}))>>;
};

template <typename Left, typename Right> struct MergeResources;
template <typename... Left>
struct MergeResources<ResourceList<Left...>, ResourceList<>> {
  using Type = ResourceList<Left...>;
};
template <typename... Left, typename First, typename... Rest>
struct MergeResources<ResourceList<Left...>, ResourceList<First, Rest...>> {
  using Type = typename MergeResources<
      typename InsertResource<ResourceList<Left...>, First>::Type,
      ResourceList<Rest...>>::Type;
};

template <typename... Lists> struct MergeAll;
template <> struct MergeAll<> {
  using Type = ResourceList<>;
};
template <typename First, typename... Rest> struct MergeAll<First, Rest...> {
  using Type =
      typename MergeResources<First, typename MergeAll<Rest...>::Type>::Type;
};

/** @brief Parameter resources of a policy; custom policies may specialize it. */
template <typename Policy> struct PolicyResources {
  using Type = ResourceList<>;
};
template <template <typename...> typename Policy, typename... Arguments>
struct PolicyResources<Policy<Arguments...>> {
  using Type =
      typename MergeAll<typename PolicyResources<Arguments>::Type...>::Type;
};
template <typename... Policies>
struct PolicyResources<std::tuple<Policies...>> {
  using Type =
      typename MergeAll<typename PolicyResources<Policies>::Type...>::Type;
};

template <typename Stage> struct StageResources {
  using Type = typename PolicyResources<typename Stage::Policies>::Type;
};
template <> struct StageResources<void> {
  using Type = ResourceList<>;
};
template <CodeEmission Emission, typename... Entries>
struct StageResources<Stage::Placed<Emission, Entries...>> {
  using Type =
      typename MergeAll<typename StageResources<Entries>::Type...>::Type;
};

template <typename Pipeline> struct PipelineResources;
template <typename Binding, typename... Entries>
struct PipelineResources<Pipeline<Binding, Entries...>> {
  using Type =
      typename MergeAll<typename StageResources<Entries>::Type...>::Type;
};

struct DiscoveryBinding {
  struct FrameState;
  using Instrumentation = NoInstrumentation;
};

template <ResourceKey Key, typename Family> struct ParameterBlock {
  Family data{};
  constexpr Family &resource() { return data; }
  constexpr const Family &resource() const { return data; }
};

template <ResourceKey Key, typename List> struct FindResource;
template <ResourceKey Key> struct FindResource<Key, ResourceList<>> {
  using Type = void;
};
template <ResourceKey Key, typename First, typename... Rest>
struct FindResource<Key, ResourceList<First, Rest...>> {
  using Type = std::conditional_t<
      First::KEY.view() == Key.view(), First,
      typename FindResource<Key, ResourceList<Rest...>>::Type>;
};
template <typename Resource> struct FamilyOf {
  using Type = typename Resource::Family;
};
template <> struct FamilyOf<void> {
  using Type = void;
};

template <typename List> struct ParameterSet;
template <typename... Resources>
struct ParameterSet<ResourceList<Resources...>>
    : ParameterBlock<Resources::KEY, typename Resources::Family>... {
  using ResourceTypes = ResourceList<Resources...>;
  template <ResourceKey Key>
  using Resource = typename FindResource<Key, ResourceTypes>::Type;
  template <ResourceKey Key>
  using Family = typename FamilyOf<Resource<Key>>::Type;
  template <ResourceKey Key>
  static constexpr bool HAS = !std::is_void_v<Resource<Key>>;

  template <ResourceKey Key>
    requires(HAS<Key>)
  constexpr auto &get() {
    return static_cast<ParameterBlock<Key, Family<Key>> &>(*this).resource();
  }
  template <ResourceKey Key>
    requires(HAS<Key>)
  constexpr const auto &get() const {
    return static_cast<const ParameterBlock<Key, Family<Key>> &>(*this)
        .resource();
  }
  template <typename Visitor>
  __attribute__((always_inline)) constexpr void visit(Visitor &&visitor) {
    (visitor.template operator()<Resources>(get<Resources::KEY>()), ...);
  }
  template <typename Visitor>
  __attribute__((always_inline)) constexpr void visit(Visitor &&visitor) const {
    (visitor.template operator()<Resources>(get<Resources::KEY>()), ...);
  }
};

template <typename Family, ResourceKind Kind> struct ClockState {};
template <typename Family> struct ClockState<Family, ResourceKind::SOURCE> {
  float primary = 0.0f;
  float secondary = 0.0f;
  float angle = 0.0f;
  float noise_time = 0.0f;
};
template <typename Family> struct RotationClock {};
template <> struct RotationClock<Warp::AffineParams> {
  float rotation = 0.0f;
};
template <typename Family>
struct ClockState<Family, ResourceKind::WARP> : RotationClock<Family> {
  float phase = 0.0f;
};
template <typename Family> struct ClockState<Family, ResourceKind::SURFACE> {
  float phase = 0.0f;
};
template <typename Resource>
struct ResourceClock : ClockState<typename Resource::Family, Resource::KIND> {};

template <typename Resource>
inline constexpr bool RESOURCE_NOISE =
    std::is_same_v<typename Resource::Family, Warp::VectorNoiseParams> ||
    std::is_same_v<typename Resource::Family,
                   Source::ProjectedNoiseSourceParams> ||
    std::is_same_v<typename Resource::Family,
                   Source::SphericalNoiseSourceParams> ||
    std::is_same_v<typename Resource::Family, Surface::SurfaceNoiseParams> ||
    std::is_same_v<typename Resource::Family, Surface::DirectSurfaceParams>;

template <typename Resource, bool = RESOURCE_NOISE<Resource>>
struct ResourceNoise {};
template <typename Resource> struct ResourceNoise<Resource, true> {
  FastNoiseLite noise;
};
template <typename Resource, bool = RESOURCE_NOISE<Resource>>
struct ResourceNoiseView {};
template <typename Resource> struct ResourceNoiseView<Resource, true> {
  const FastNoiseLite *noise = nullptr;
};
template <typename Resource>
struct ResourceFrame : ResourceClock<Resource>, ResourceNoiseView<Resource> {};

template <template <typename> typename Storage, typename Resources>
struct ResourceStorage;
template <template <typename> typename Storage, typename... Resources>
struct ResourceStorage<Storage, ResourceList<Resources...>>
    : Storage<Resources>... {
  using ResourceTypes = ResourceList<Resources...>;
  template <ResourceKey Key> auto &get() {
    using Resource = typename FindResource<Key, ResourceTypes>::Type;
    static_assert(!std::is_void_v<Resource>,
                  "composed runtime resource does not exist");
    return static_cast<Storage<Resource> &>(*this);
  }
  template <ResourceKey Key> const auto &get() const {
    using Resource = typename FindResource<Key, ResourceTypes>::Type;
    static_assert(!std::is_void_v<Resource>,
                  "composed runtime resource does not exist");
    return static_cast<const Storage<Resource> &>(*this);
  }
};

} // namespace ComposedDetail

/** @brief Present-only parameter storage derived from a declared ranked pipeline. */
template <typename Spec>
using ParamsFor = ComposedDetail::ParameterSet<
    typename ComposedDetail::PipelineResources<typename Spec::template Pipeline<
        ComposedDetail::DiscoveryBinding>>::Type>;

} // namespace Pullback
