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

#include "core/math/noise_field.h"
#include "core/render/pullback/lens.h"
#include "core/render/pullback/source.h"
#include "core/render/pullback/sphere_animation.h"
#include "core/render/pullback/stage.h"
#include "core/render/pullback/warp.h"

namespace Pullback {

/** @brief Compile-time identity of one composed parameter/resource instance. */
template <size_t N> struct ResourceKey {
  char text[N]; ///< Key characters, including the terminating NUL.
  /**
   * @brief Captures a string literal as a key.
   * @param value NUL-terminated key literal.
   */
  constexpr ResourceKey(const char (&value)[N]) {
    for (size_t i = 0; i < N; ++i)
      text[i] = value[i];
  }
  /** @brief The key text without its terminating NUL.
   *  @return View over `text`. */
  constexpr std::string_view view() const { return {text, N - 1}; }
};

/** @brief Pipeline role of a composed parameter resource. */
enum class ResourceKind {
  SOURCE,     ///< Scalar pattern source.
  PROJECTION, ///< Sphere-to-plane projection.
  WARP,       ///< Planar warp.
  SURFACE,    ///< Displacement surface.
  LENS,       ///< Mobius lens.
  VALUE,      ///< Value transfer and coverage.
  COLOR,      ///< Colorizer.
  MODULATION  ///< Ordered color modulation.
};

/** @brief The storage required by one provider instance. */
template <ResourceKey Key, typename FamilyT, ResourceKind KindV>
struct ParameterResource {
  static constexpr auto KEY = Key;            ///< Instance key.
  static constexpr ResourceKind KIND = KindV; ///< Pipeline role.
  using Family = FamilyT; ///< Parameter family stored for this instance.
  /// Standard instance keys, in canonical storage order.
  static constexpr std::string_view INSTANCE_ORDER[] = {
      "source", "projection", "outer_warp", "inner_warp",   "surface",
      "lens",   "value",      "color",      "displacement", "modulation"};
  /// Index of `Key` in `INSTANCE_ORDER`, or its size for a non-standard
  /// key.
  static constexpr size_t ORDER = [] {
    for (size_t i = 0; i < std::size(INSTANCE_ORDER); ++i)
      if (Key.view() == INSTANCE_ORDER[i])
        return i;
    return std::size(INSTANCE_ORDER);
  }();
  /// Whether `Key` is one of `INSTANCE_ORDER`.
  static constexpr bool STANDARD = ORDER < std::size(INSTANCE_ORDER);
};

namespace ComposedDetail {

/** @brief Type list of ParameterResource types, ordered by `ORDER`. */
template <typename... T> struct ResourceList {};

/** @brief Inserts `Resource` into `List` in `ORDER`, dropping a same-key
 *  duplicate. */
template <typename List, typename Resource> struct InsertResource;
/** @brief Inserting into an empty list yields a one-element list. */
template <typename Resource> struct InsertResource<ResourceList<>, Resource> {
  using Type = ResourceList<Resource>;
};
/** @brief Recursive case of InsertResource. */
template <typename First, typename... Rest, typename Resource>
struct InsertResource<ResourceList<First, Rest...>, Resource> {
  /// Whether `First` and `Resource` share a key.
  static constexpr bool SAME_KEY = First::KEY.view() == Resource::KEY.view();
  static_assert(!SAME_KEY || std::is_same_v<First, Resource>,
                "composed resource key has conflicting families or roles");
  /** @brief Type-level prepend of `First`; declared only for `decltype`.
   *  @tparam T Elements of the list to prepend to.
   *  @return The list with `First` in front. */
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

/** @brief Sorted union of two resource lists. */
template <typename Left, typename Right> struct MergeResources;
/** @brief Merging an empty list yields `Left`. */
template <typename... Left>
struct MergeResources<ResourceList<Left...>, ResourceList<>> {
  using Type = ResourceList<Left...>;
};
/** @brief Recursive case of MergeResources. */
template <typename... Left, typename First, typename... Rest>
struct MergeResources<ResourceList<Left...>, ResourceList<First, Rest...>> {
  using Type = typename MergeResources<
      typename InsertResource<ResourceList<Left...>, First>::Type,
      ResourceList<Rest...>>::Type;
};

/** @brief Sorted union of any number of resource lists. */
template <typename... Lists> struct MergeAll;
/** @brief The union of no lists is empty. */
template <> struct MergeAll<> {
  using Type = ResourceList<>;
};
/** @brief Recursive case of MergeAll. */
template <typename First, typename... Rest> struct MergeAll<First, Rest...> {
  using Type =
      typename MergeResources<First, typename MergeAll<Rest...>::Type>::Type;
};

/** @brief Parameter resources of a policy; custom policies may specialize it. */
template <typename Policy> struct PolicyResources {
  using Type = ResourceList<>; ///< No resources by default.
};
/** @brief A policy template's resources: the union over its arguments. */
template <template <typename...> typename Policy, typename... Arguments>
struct PolicyResources<Policy<Arguments...>> {
  using Type =
      typename MergeAll<typename PolicyResources<Arguments>::Type...>::Type;
};
/** @brief A policy tuple's resources: the union over its policies. */
template <typename... Policies>
struct PolicyResources<std::tuple<Policies...>> {
  using Type =
      typename MergeAll<typename PolicyResources<Policies>::Type...>::Type;
};

/** @brief Parameter resources of a stage's `Policies`. */
template <typename Stage> struct StageResources {
  /// Union of the stage's policy resources.
  using Type = typename PolicyResources<typename Stage::Policies>::Type;
};
/** @brief An absent stage needs no resources. */
template <> struct StageResources<void> {
  using Type = ResourceList<>;
};
/** @brief A placed stage group's resources: the union over its entries. */
template <CodeEmission Emission, typename... Entries>
struct StageResources<Stage::Placed<Emission, Entries...>> {
  using Type =
      typename MergeAll<typename StageResources<Entries>::Type...>::Type;
};

/** @brief Parameter resources of a whole pipeline. */
template <typename Pipeline> struct PipelineResources;
/** @brief A pipeline's resources: the union over its stages. */
template <typename Binding, typename... Entries>
struct PipelineResources<Pipeline<Binding, Entries...>> {
  using Type =
      typename MergeAll<typename StageResources<Entries>::Type...>::Type;
};

/** @brief Declaration binding for resource and topology discovery. */
struct DiscoveryBinding {
  struct FrameState;
  using Instrumentation = NoInstrumentation; ///< Uninstrumented.
};

/** @brief One keyed parameter family, as a ParameterSet base.
 * @tparam Key Resource key.
 * @tparam Family Parameter record. */
template <ResourceKey Key, typename Family> struct ParameterBlock {
  Family data{}; ///< The family's parameter values.
  /** @brief Mutable access.
   *  @return `data`. */
  constexpr Family &resource() { return data; }
  /** @brief Const access.
   *  @return `data`. */
  constexpr const Family &resource() const { return data; }
};

/** @brief The resource in `List` with key `Key`, or void. */
template <ResourceKey Key, typename List> struct FindResource;
/** @brief No match in an empty list. */
template <ResourceKey Key> struct FindResource<Key, ResourceList<>> {
  using Type = void;
};
/** @brief Recursive case of FindResource. */
template <ResourceKey Key, typename First, typename... Rest>
struct FindResource<Key, ResourceList<First, Rest...>> {
  using Type = std::conditional_t<
      First::KEY.view() == Key.view(), First,
      typename FindResource<Key, ResourceList<Rest...>>::Type>;
};
/** @brief The parameter family of a resource, or void for void. */
template <typename Resource> struct FamilyOf {
  /// `Resource::Family`.
  using Type = typename Resource::Family;
};
/** @brief A missing resource has no family. */
template <> struct FamilyOf<void> {
  using Type = void;
};

template <typename List> struct ParameterSet;
/** @brief Keyed parameter storage holding one family per resource.
 * @tparam Resources Parameter resources. */
template <typename... Resources>
struct ParameterSet<ResourceList<Resources...>>
    : ParameterBlock<Resources::KEY, typename Resources::Family>... {
  using ResourceTypes = ResourceList<Resources...>; ///< The stored resources.
  /// The resource keyed `Key`, or void.
  template <ResourceKey Key>
  using Resource = typename FindResource<Key, ResourceTypes>::Type;
  /// The family keyed `Key`, or void.
  template <ResourceKey Key>
  using Family = typename FamilyOf<Resource<Key>>::Type;
  /// Whether a resource keyed `Key` is stored.
  template <ResourceKey Key>
  static constexpr bool HAS = !std::is_void_v<Resource<Key>>;

  /** @brief The family keyed `Key`.
   *  @tparam Key Resource key.
   *  @return Mutable reference to the family. */
  template <ResourceKey Key>
    requires(HAS<Key>)
  constexpr auto &get() {
    return static_cast<ParameterBlock<Key, Family<Key>> &>(*this).resource();
  }
  /** @brief The family keyed `Key`.
   *  @tparam Key Resource key.
   *  @return Const reference to the family. */
  template <ResourceKey Key>
    requires(HAS<Key>)
  constexpr const auto &get() const {
    return static_cast<const ParameterBlock<Key, Family<Key>> &>(*this)
        .resource();
  }
  /** @brief Calls `visitor.operator()<Resource>(family)` for each resource
   *  in storage order.
   *  @tparam Visitor Callable with a templated call operator.
   *  @param visitor Receives each mutable family. */
  template <typename Visitor>
  __attribute__((always_inline)) constexpr void visit(Visitor &&visitor) {
    (visitor.template operator()<Resources>(get<Resources::KEY>()), ...);
  }
  /** @brief Calls `visitor.operator()<Resource>(family)` for each resource
   *  in storage order.
   *  @tparam Visitor Callable with a templated call operator.
   *  @param visitor Receives each const family. */
  template <typename Visitor>
  __attribute__((always_inline)) constexpr void visit(Visitor &&visitor) const {
    (visitor.template operator()<Resources>(get<Resources::KEY>()), ...);
  }
};

/** @brief Per-frame clocks of one resource; empty for unclocked kinds. */
template <typename Family, ResourceKind Kind> struct ClockState {};
/** @brief Source pattern clocks.
 * @details Values are equivalent modulo their periods; column-time samples
 * may lie outside the canonical ranges. */
template <typename Family> struct ClockState<Family, ResourceKind::SOURCE> {
  float primary = 0.0f;    ///< Primary phase, radians modulo 2pi.
  float secondary = 0.0f;  ///< Secondary phase, radians modulo 2pi.
  float angle = 0.0f;      ///< Pattern rotation, radians modulo 2pi.
  float noise_time = 0.0f; ///< Noise time coordinate modulo 1.
};
/** @brief Extra rotation clock of a warp family; empty by default. */
template <typename Family> struct RotationClock {};
/** @brief The affine warp's accumulated rotation.
 * @details Column-time samples may lie outside [0, 2pi). */
template <> struct RotationClock<Warp::AffineParams> {
  float rotation = 0.0f; ///< Frame rotation, radians modulo 2pi.
};
/** @brief Warp phase clock.
 * @details Column-time samples may lie outside [0, 1). */
template <typename Family>
struct ClockState<Family, ResourceKind::WARP> : RotationClock<Family> {
  float phase = 0.0f; ///< Warp phase modulo 1.
};
/** @brief Surface phase clock.
 * @details Values are equivalent modulo their periods; non-ripple column-time
 * samples may lie outside [0, 1). */
template <typename Family> struct ClockState<Family, ResourceKind::SURFACE> {
  float phase = 0.0f; ///< Phase modulo 1; ripple frame count in [0, period).
};
/** @brief Phase of the animation generating a Mobius lens's b. */
template <> struct ClockState<Lens::MobiusLensParams, ResourceKind::LENS> {
  float angle = 0.0f;     ///< Generator angle, radians.
  float magnitude = 0.0f; ///< Generator magnitude; zero when unanimated.
  bool generated = false; ///< Whether the frame's b is the generator's.
};
/** @brief Ripple burst clocks and active wave ages. */
template <>
struct ClockState<SphereAnimation::BurstRippleParams, ResourceKind::SURFACE>
    : SphereAnimation::RippleState {};
/** @brief Modulation source clock.
 * @tparam Family Authored modulation family. */
template <typename Family> struct ClockState<Family, ResourceKind::MODULATION> {
  Modulation::State modulation_clock; ///< Selected source clock.
};
/** @brief Spatial source storage of an authored layer.
 * @tparam Resource Keyed resource.
 * @tparam Spatial Whether the layer samples a spatial source. */
template <typename Resource,
          bool Spatial = Resource::KIND == ResourceKind::MODULATION>
struct ResourceModulation {};
/** @brief Modulation resource storage.
 * @tparam Resource Authored modulation resource. */
template <typename Resource> struct ResourceModulation<Resource, true> {
  static constexpr bool SPATIAL =
      Resource::Family::SOURCE !=
      Pattern::SourceKind::DISPLACEMENT; ///< Spatial capability.
  std::conditional_t<SPATIAL,
                     Modulation::Resources<Resource::Family::EVALUATION ==
                                           Pattern::Evaluation::BAKED_CUBE>,
                     std::monostate>
      modulation_resources; ///< Source cache.
};
/** @brief Borrowed modulation endpoint pair.
 * @tparam Resource Keyed resource.
 * @tparam Enabled Whether this resource modulates color. */
template <typename Resource,
          bool Enabled = Resource::KIND == ResourceKind::MODULATION>
struct ResourceModulationView {};
/** @brief Modulation endpoint pair.
 * @tparam Resource Keyed modulation resource. */
template <typename Resource> struct ResourceModulationView<Resource, true> {
  const Modulation::Prepared *modulation = nullptr; ///< Current source.
  const Modulation::Prepared *previous_modulation =
      nullptr; ///< Retained source.
};
/** @brief Retained spatial source storage.
 * @tparam Resource Keyed resource.
 * @tparam Enabled Whether this resource owns a source. */
template <typename Resource,
          bool Enabled = Resource::KIND == ResourceKind::MODULATION>
struct ResourceModulationHistory {};
/** @brief Retention specialized for a modulation resource.
 * @tparam Resource Keyed modulation resource. */
template <typename Resource> struct ResourceModulationHistory<Resource, true> {
  static constexpr bool SPATIAL =
      Resource::Family::SOURCE !=
      Pattern::SourceKind::DISPLACEMENT; ///< Source capability.
  std::conditional_t<SPATIAL, Modulation::Prepared, std::monostate>
      modulation; ///< Retained inputs.
  std::array<int8_t, SPATIAL && Resource::Family::EVALUATION ==
                                    Pattern::Evaluation::BAKED_CUBE
                         ? Pattern::SphereLutView::SIZE
                         : 0>
      table{}; ///< Retained cube.
};

/** @brief The clocks of `Resource`. */
template <typename Resource>
struct ResourceClock : ClockState<typename Resource::Family, Resource::KIND> {};

/// Whether `Resource`'s family samples a FastNoiseLite field.
template <typename Resource>
inline constexpr bool RESOURCE_NOISE =
    std::is_same_v<typename Resource::Family, Warp::VectorNoiseParams> ||
    std::is_same_v<typename Resource::Family,
                   Source::ProjectedNoiseSourceParams> ||
    std::is_same_v<typename Resource::Family,
                   Source::SphericalNoiseSourceParams> ||
    std::is_same_v<typename Resource::Family, Surface::SurfaceNoiseParams> ||
    std::is_same_v<typename Resource::Family, Surface::DirectSurfaceParams>;

/** @brief Persistent ripple scheduler, empty for other families.
 * @tparam Resource Keyed resource.
 * @tparam Enabled Whether the resource schedules ripple bursts. */
template <typename Resource,
          bool Enabled = std::is_same_v<typename Resource::Family,
                                        SphereAnimation::BurstRippleParams>>
struct ResourceRippleScheduler {};
/** @brief Persistent ripple burst scheduler.
 * @tparam Resource Ripple resource. */
template <typename Resource> struct ResourceRippleScheduler<Resource, true> {
  SphereAnimation::RippleState ripples; ///< Active bursts and clock.
};
/** @brief Member clocks, excluding arena-owned burst schedulers.
 * @tparam Resource Keyed resource.
 * @tparam Arena Whether the clock is arena-owned. */
template <typename Resource,
          bool Arena = std::is_same_v<typename Resource::Family,
                                      SphereAnimation::BurstRippleParams>>
struct ResourceRuntimeClock : ResourceClock<Resource> {};
/** @brief Empty member slot of an arena-owned burst scheduler.
 * @tparam Resource Ripple resource. */
template <typename Resource> struct ResourceRuntimeClock<Resource, true> {};
/** @brief Owned resources, including noise where required. */
template <typename Resource, bool = RESOURCE_NOISE<Resource>>
struct ResourceNoise : ResourceModulation<Resource>,
                       ResourceRippleScheduler<Resource> {};
/** @brief Noise resource: owns its generator. */
template <typename Resource>
struct ResourceNoise<Resource, true> : ResourceModulation<Resource> {
  FastNoiseLite noise; ///< Seeded generator.
};
/** @brief Non-owning view of a resource's noise; empty otherwise. */
template <typename Resource, bool = RESOURCE_NOISE<Resource>>
struct ResourceNoiseView {};
/** @brief Noise resource: points at the runtime's generator. */
template <typename Resource> struct ResourceNoiseView<Resource, true> {
  const FastNoiseLite *noise = nullptr; ///< Runtime-owned generator.
};
/** @brief Optional persistent radial table for one resource.
 * @tparam Resource Keyed resource.
 * @tparam Enabled Whether the effect bakes this resource's ripple curves. */
template <typename Resource, bool Enabled> struct ResourceRippleRadialTable {};
/** @brief Enabled persistent ripple curves.
 * @tparam Resource Burst-ripple resource. */
template <typename Resource> struct ResourceRippleRadialTable<Resource, true> {
  SphereAnimation::RadialTables radial_tables; ///< Frame endpoint curves.
};

/** @brief Optional borrowed ripple transfer curves.
 * @tparam Resource Keyed resource.
 * @tparam Enabled Whether this family schedules burst ripples. */
template <typename Resource,
          bool Enabled = std::is_same_v<typename Resource::Family,
                                        SphereAnimation::BurstRippleParams>>
struct ResourceRippleRadialView {};
/** @brief Borrowed ripple curves for a column time.
 * @tparam Resource Burst-ripple resource. */
template <typename Resource> struct ResourceRippleRadialView<Resource, true> {
  const SphereAnimation::RadialTables *radial_tables =
      nullptr; ///< Runtime-owned curves.
  float radial_fraction =
      -1.0f; ///< Current endpoint weight; negative disables.
};

/** @brief Frame snapshot of one resource: its clocks and noise view. */
template <typename Resource>
struct ResourceFrame : ResourceClock<Resource>,
                       ResourceNoiseView<Resource>,
                       ResourceRippleRadialView<Resource>,
                       ResourceModulationView<Resource> {};

/** @brief One `Storage<Resource>` per resource, addressed by key. */
template <template <typename> typename Storage, typename Resources>
struct ResourceStorage;
/** @brief ResourceStorage over a resource list. */
template <template <typename> typename Storage, typename... Resources>
struct ResourceStorage<Storage, ResourceList<Resources...>>
    : Storage<Resources>... {
  using ResourceTypes = ResourceList<Resources...>; ///< The stored resources.
  /** @brief Storage of the resource keyed `Key`.
   *  @tparam Key Resource key; must exist.
   *  @return Mutable reference to its `Storage`. */
  template <ResourceKey Key> auto &get() {
    using Resource = typename FindResource<Key, ResourceTypes>::Type;
    static_assert(!std::is_void_v<Resource>,
                  "composed runtime resource does not exist");
    return static_cast<Storage<Resource> &>(*this);
  }
  /** @brief Storage of the resource keyed `Key`.
   *  @tparam Key Resource key; must exist.
   *  @return Const reference to its `Storage`. */
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
