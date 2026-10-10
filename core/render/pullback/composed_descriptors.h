/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/** @file composed_descriptors.h
 * @brief Compile-time stage classification, Spec checks and parameter
 * counting for composed effects.
 */

namespace ComposedDetail {

/// Slider names of the Mobius coefficients, in `math::MobiusParams` order.
inline constexpr const char *MOBIUS_PARAM_NAMES[] = {
    "Mobius A Re", "Mobius A Im", "Mobius B Re", "Mobius B Im",
    "Mobius C Re", "Mobius C Im", "Mobius D Re", "Mobius D Im"};
/**
 * @brief Slider name of a warp instance's speed parameter.
 * @param key Warp instance key.
 * @return The numbered name for `outer_warp` / `inner_warp`, else the plain
 *         name.
 */
constexpr std::string_view warp_speed_name(std::string_view key) {
  return key == "outer_warp"   ? "Planar Warp 1 Speed"
         : key == "inner_warp" ? "Planar Warp 2 Speed"
                               : "Planar Warp Speed";
}

/** @brief Whether @p T is a `Stage::Sample` stage. @tparam T Stage type. */
template <typename T> struct IsSampleStage : std::false_type {};
template <typename S, typename W, typename C>
struct IsSampleStage<Stage::Sample<S, W, C>> : std::true_type {};
/**
 * @brief The `ProjectionCoverageMode` a `ProjectionCoverage` policy implements.
 * @tparam T Coverage policy type.
 */
template <typename T> struct ProjectionCoverageModeOf {
  /** Out-of-range sentinel: @p T is not a known coverage policy. */
  static constexpr auto VALUE = static_cast<ProjectionCoverageMode>(255);
};
/** @brief `ProjectionCoverage::None` is `NONE`. */
template <> struct ProjectionCoverageModeOf<ProjectionCoverage::None> {
  static constexpr auto VALUE = ProjectionCoverageMode::NONE; ///< The mode.
};
/** @brief `ProjectionCoverage::Weight` is `WEIGHT`. */
template <> struct ProjectionCoverageModeOf<ProjectionCoverage::Weight> {
  static constexpr auto VALUE = ProjectionCoverageMode::WEIGHT; ///< The mode.
};
/** @brief `ProjectionCoverage::WeightSquared` is `WEIGHT_SQUARED`. */
template <> struct ProjectionCoverageModeOf<ProjectionCoverage::WeightSquared> {
  static constexpr auto VALUE =
      ProjectionCoverageMode::WEIGHT_SQUARED; ///< The mode.
};
/**
 * @brief `ProjectionCoverage::EdgeFade` is `EDGE_FADE`.
 * @tparam P Edge-width provider.
 */
template <typename P>
struct ProjectionCoverageModeOf<ProjectionCoverage::EdgeFade<P>> {
  static constexpr auto VALUE =
      ProjectionCoverageMode::EDGE_FADE; ///< The mode.
};

/** @brief Whether @p T is a `Stage::Lens` stage. @tparam T Stage type. */
template <typename T> struct IsLensStage : std::false_type {};
template <typename P> struct IsLensStage<Stage::Lens<P>> : std::true_type {};
/** @brief Whether @p T is a `Stage::Displace` stage. @tparam T Stage type. */
template <typename T> struct IsSurfaceStage : std::false_type {};
template <typename P>
struct IsSurfaceStage<Stage::Displace<P>> : std::true_type {};
/** @brief Whether @p T is a `Stage::Project` stage. @tparam T Stage type. */
template <typename T> struct IsProjectStage : std::false_type {};
template <typename P>
struct IsProjectStage<Stage::Project<P>> : std::true_type {};
/** @brief Whether @p T is a `Stage::Transfer` stage. @tparam T Stage type. */
template <typename T> struct IsTransferStage : std::false_type {};
template <typename P>
struct IsTransferStage<Stage::Transfer<P>> : std::true_type {};
/** @brief Whether @p T is a `Stage::ApplyCoverage` stage. @tparam T Stage type. */
template <typename T> struct IsCoverageStage : std::false_type {};
template <typename P>
struct IsCoverageStage<Stage::ApplyCoverage<P>> : std::true_type {};
/** @brief Whether @p T is a `Lens::Mobius` policy. @tparam T Policy type. */
template <typename T> struct IsMobiusLens : std::false_type {};
template <typename P> struct IsMobiusLens<Lens::Mobius<P>> : std::true_type {};

/**
 * @brief Index of the first stage of @p Pipeline matching @p Predicate.
 * @tparam Pipeline Stage pipeline.
 * @tparam Predicate Stage trait with a boolean `value`.
 * @tparam Index Stage the search starts at.
 * @return The stage index, or `Pipeline::STAGE_COUNT` when none matches.
 */
template <typename Pipeline, template <typename> class Predicate,
          size_t Index = 0>
consteval size_t stage_index() {
  if constexpr (Index == Pipeline::STAGE_COUNT)
    return Index;
  else if constexpr (Predicate<
                         typename Pipeline::template stage_at<Index>>::value)
    return Index;
  else
    return stage_index<Pipeline, Predicate, Index + 1>();
}

/**
 * @brief Whether @p Policy, or a provider it wraps, tracks path length.
 * @tparam Policy Stage policy.
 */
template <typename Policy> struct PathTracked : std::false_type {};
template <template <typename...> class Policy, typename... Arguments>
struct PathTracked<Policy<Arguments...>>
    : std::disjunction<PathTracked<Arguments>...> {};
template <typename B, ResourceKey Key, typename Family, bool Track,
          ResourceKey SourceKey>
struct PathTracked<WarpProvider<B, Key, Family, Track, SourceKey>>
    : std::bool_constant<Track> {};
template <typename B, typename Family, bool Track, ResourceKey Key>
struct PathTracked<SurfaceProvider<B, Family, Track, Key>>
    : std::bool_constant<Track> {};
template <typename Provider, typename Mode, uint8_t Harmonic>
struct PathTracked<Warp::PolarChart<Provider, Mode, Harmonic>>
    : PathTracked<Provider> {};
template <typename Provider, math::NoiseBasis Basis, typename Envelope>
struct PathTracked<Warp::VectorNoise<Provider, Basis, Envelope>>
    : PathTracked<Provider> {};
template <typename Provider, math::NoiseBasis Basis, typename Integrator,
          math::TangentLimit Limit>
struct PathTracked<Surface::CurlNoise<Provider, Basis, Integrator, Limit>>
    : PathTracked<Provider> {};
template <typename Provider, math::NoiseBasis Basis>
struct PathTracked<Surface::DirectNoise<Provider, Basis>>
    : PathTracked<Provider> {};
/** @brief Whether any policy of @p Stage tracks path length.
 * @tparam Stage Pipeline stage. */
template <typename Stage>
struct StagePathTracked : PathTracked<typename Stage::Policies> {};

/**
 * @brief Checks of a Spec's metadata fields against its pipeline.
 * @tparam Spec The effect's `Pullback::Spec`.
 * @tparam Binding Frame binding of the effect.
 * @tparam Pipeline The pipeline the Spec declares for @p Binding.
 */
template <typename Spec, typename Binding, typename Pipeline>
struct PipelineMetadata {
  /** Whether any stage's warp or surface provider tracks path length. */
  static constexpr bool PATH_TRACKED =
      Pipeline::template any_stage<StagePathTracked>;
  /// The pipeline's `Stage::Lens`, or void when it has none.
  using LensStage = typename Pipeline::template stage_matching<IsLensStage>;
  /// The pipeline's `Stage::Project`, or void when it has none.
  using ProjectStage =
      typename Pipeline::template stage_matching<IsProjectStage>;
  /// The pipeline's `Stage::Sample`, or void when it has none.
  using SampleStage = typename Pipeline::template stage_matching<IsSampleStage>;
  /** Whether `Spec::COVERAGE` names the Sample stage's coverage policy
      (`NONE` without a Sample stage). */
  static constexpr bool COVERAGE_MATCHES = [] {
    if constexpr (std::is_void_v<SampleStage>)
      return Spec::COVERAGE == ProjectionCoverageMode::NONE;
    else
      return ProjectionCoverageModeOf<
                 typename SampleStage::CoveragePolicy>::VALUE == Spec::COVERAGE;
  }();
  /** Whether `Spec::LensPolicy` is the lens stage's policy, or void for no
      lens or a Mobius lens. */
  static constexpr bool LENS_MATCHES = [] {
    if constexpr (std::is_void_v<LensStage>)
      return std::is_void_v<typename Spec::LensPolicy>;
    else if constexpr (IsMobiusLens<typename LensStage::LensPolicy>::value)
      return std::is_void_v<typename Spec::LensPolicy>;
    else
      return std::is_same_v<typename LensStage::LensPolicy,
                            typename Spec::LensPolicy>;
  }();
  /** Whether the Project stage uses the `Spec::PROJECTION` policy; false
      without a Project stage. */
  static constexpr bool PROJECTION_MATCHES = [] {
    if constexpr (std::is_void_v<ProjectStage>)
      return false;
    else
      return std::is_same_v<
          typename ProjectStage::ProjectionPolicy,
          typename ProjectionPolicyFor<Spec::PROJECTION, Binding>::Type>;
  }();
  /** Whether the surface and lens stage order agrees with
      `Spec::SURFACE_PLACEMENT`; true unless both stages are present. */
  static constexpr bool SURFACE_PLACEMENT_MATCHES = [] {
    constexpr size_t SURFACE = stage_index<Pipeline, IsSurfaceStage>();
    constexpr size_t LENS = stage_index<Pipeline, IsLensStage>();
    if constexpr (SURFACE == Pipeline::STAGE_COUNT ||
                  LENS == Pipeline::STAGE_COUNT)
      return true;
    else
      return (SURFACE < LENS) ==
             (Spec::SURFACE_PLACEMENT == SurfacePlacement::BEFORE_LENS);
  }();
};

/** @brief The outer camera reads the "projection" family. @tparam B Frame
    binding. */
template <typename B> struct PolicyResources<OuterCameraProvider<B>> {
  using Type = ResourceList<ParameterResource<"projection", ProjectionParams,
                                              ResourceKind::PROJECTION>>;
};
template <typename B>
struct PolicyResources<ProjectionProvider<B>>
    : PolicyResources<OuterCameraProvider<B>> {};
/**
 * @brief A lens provider declares a Mobius lens family under @p Key.
 * @tparam B Frame binding.
 * @tparam Key Parameter resource key.
 */
template <typename B, ResourceKey Key>
struct PolicyResources<LensProvider<B, Key>> {
  using Type = ResourceList<
      ParameterResource<Key, MobiusLensParams, ResourceKind::LENS>>;
};
/**
 * @brief A warp provider declares its warp family under @p Key.
 * @tparam B Frame binding.
 * @tparam Key Parameter resource key.
 * @tparam Family Warp parameter family.
 * @tparam Track Whether path length is tracked.
 * @tparam SourceKey Resource key of the source the warp reads.
 */
template <typename B, ResourceKey Key, typename Family, bool Track,
          ResourceKey SourceKey>
struct PolicyResources<WarpProvider<B, Key, Family, Track, SourceKey>> {
  using Type = ResourceList<ParameterResource<Key, Family, ResourceKind::WARP>>;
};
/**
 * @brief A surface provider declares its surface family under @p Key.
 * @tparam B Frame binding.
 * @tparam Family Surface parameter family.
 * @tparam Track Whether path length is tracked.
 * @tparam Key Parameter resource key.
 */
template <typename B, typename Family, bool Track, ResourceKey Key>
struct PolicyResources<SurfaceProvider<B, Family, Track, Key>> {
  using Type =
      ResourceList<ParameterResource<Key, Family, ResourceKind::SURFACE>>;
};
/**
 * @brief A source provider declares its source family under @p Key.
 * @tparam B Frame binding.
 * @tparam Family Source parameter family.
 * @tparam Key Parameter resource key.
 */
template <typename B, typename Family, ResourceKey Key>
struct PolicyResources<SourceProvider<B, Family, Key>> {
  using Type =
      ResourceList<ParameterResource<Key, Family, ResourceKind::SOURCE>>;
};
/**
 * @brief A value provider declares its value family under @p Key.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 * @tparam Key Parameter resource key.
 */
template <typename B, typename Family, ResourceKey Key>
struct PolicyResources<ValueProvider<B, Family, Key>> {
  using Type =
      ResourceList<ParameterResource<Key, Family, ResourceKind::VALUE>>;
};
/**
 * @brief A colour provider declares the "color" family.
 * @tparam B Frame binding.
 * @tparam Hue Hue rotation mode.
 * @tparam Brightness Brightness envelope.
 */
template <typename B, HueMode Hue, Color::BrightnessEnvelope Brightness>
struct PolicyResources<ColorProvider<B, Hue, Brightness>> {
  using Type = ResourceList<
      ParameterResource<"color", ColorParams, ResourceKind::COLOR>>;
};

template <typename Provider, Projection::GnomonicHemisphere Hemisphere>
struct PolicyResources<Projection::Gnomonic<Provider, Hemisphere>>
    : PolicyResources<Provider> {};
template <typename Provider, math::NoiseBasis Basis>
struct PolicyResources<Source::ProjectedNoise<Provider, Basis>>
    : PolicyResources<Provider> {};
template <typename Provider, math::NoiseBasis Basis>
struct PolicyResources<Source::SphericalNoise<Provider, Basis>>
    : PolicyResources<Provider> {};
template <typename Provider, typename Mode, uint8_t Harmonic>
struct PolicyResources<Warp::PolarChart<Provider, Mode, Harmonic>>
    : PolicyResources<Provider> {};
template <typename Provider, math::NoiseBasis Basis, typename Envelope>
struct PolicyResources<Warp::VectorNoise<Provider, Basis, Envelope>>
    : PolicyResources<Provider> {};
template <typename Provider, math::NoiseBasis Basis, typename Integrator,
          math::TangentLimit Limit>
struct PolicyResources<Surface::CurlNoise<Provider, Basis, Integrator, Limit>>
    : PolicyResources<Provider> {};
template <typename Provider, math::NoiseBasis Basis>
struct PolicyResources<Surface::DirectNoise<Provider, Basis>>
    : PolicyResources<Provider> {};

/**
 * @brief Whether a field gate is live under the Spec's projection metadata.
 * @tparam Spec Composed-effect Spec.
 * @param gate Field gate.
 * @return True when the gated field gets a slider.
 */
template <typename Spec> constexpr bool field_gate_open(FieldGate gate) {
  switch (gate) {
  case FieldGate::ALWAYS:
    return true;
  case FieldGate::ANIMATED_PROJECTION:
    return Spec::ANIMATED_PROJECTION;
  case FieldGate::CENTRAL_MERIDIAN:
    return uses_central_meridian(Spec::PROJECTION);
  case FieldGate::SINGULARITY_FADE:
    return uses_singularity_fade(Spec::PROJECTION);
  }
  return false;
}

/** Declared only: a call makes the enclosing constant evaluation ill-formed. */
void color_gate_names_unknown_topology();

/**
 * @brief Whether a colour field's topology gate is live under the Spec's
 *        `HUE` and `BRIGHTNESS` selections.
 * @tparam Spec Composed-effect Spec.
 * @param gate Topology gate; a null `field` is always open.
 * @return True when the gated colour field gets a slider.
 */
template <typename Spec> consteval bool color_gate_open(TopologyGate gate) {
  if (gate.field == nullptr)
    return true;
  const std::string_view topology = gate.field;
  unsigned value = 0;
  if (topology == Color::HUE_ROTATION_GATE.field)
    value = static_cast<unsigned>(Spec::HUE);
  else if (topology == Color::BRIGHTNESS_ENVELOPE_GATE.field)
    value = static_cast<unsigned>(Spec::BRIGHTNESS);
  else
    color_gate_names_unknown_topology();
  return (gate.values >> value) & 1U;
}

/**
 * @brief Number of sliders one resource contributes.
 * @tparam Spec Composed-effect Spec, for gate evaluation.
 * @tparam Resource Parameter resource.
 * @return Open fields, plus the speed slider of a warp or the palette-mapping
 *         slider of a colour family; the Mobius coefficients for a lens.
 */
template <typename Spec, typename Resource>
consteval size_t resource_parameter_count() {
  using Family = typename Resource::Family;
  if constexpr (Resource::KIND == ResourceKind::LENS)
    return std::size(MOBIUS_PARAM_NAMES);
  else {
    size_t count = Resource::KIND == ResourceKind::WARP ||
                           Resource::KIND == ResourceKind::COLOR
                       ? 1
                       : 0;
    for (const auto &field : Family::FIELDS) {
      if constexpr (Resource::KIND == ResourceKind::COLOR) {
        if (color_gate_open<Spec>(field.topology_gate))
          ++count;
      } else if (field.name != nullptr && field_gate_open<Spec>(field.gate))
        ++count;
    }
    return count;
  }
}
/**
 * @brief Total slider count of a resource list.
 * @tparam Spec Composed-effect Spec.
 * @tparam Resources Parameter resources.
 * @return Sum of `resource_parameter_count` over @p Resources.
 */
template <typename Spec, typename... Resources>
consteval size_t parameter_count(ResourceList<Resources...>) {
  return (resource_parameter_count<Spec, Resources>() + ... + 0);
}

/**
 * @brief Number of resources in a list whose family is @p Family.
 * @tparam Family Parameter family type.
 * @tparam Resources Parameter resources.
 * @return The instance count.
 */
template <typename Family, typename... Resources>
consteval size_t family_instances(ResourceList<Resources...>) {
  return (size_t(std::is_same_v<Family, typename Resources::Family>) + ... + 0);
}

/**
 * @brief Whether two distinct same-kind resources share a field name.
 * @tparam Resource Resource checked.
 * @tparam Other Resource compared against.
 * @return False for identical resources, different kinds and lenses.
 */
template <typename Resource, typename Other>
consteval bool resource_names_overlap() {
  if constexpr (Resource::KIND != Other::KIND ||
                std::is_same_v<Resource, Other> ||
                Resource::KIND == ResourceKind::LENS)
    return false;
  else {
    for (const auto &field : Resource::Family::FIELDS)
      for (const auto &other : Other::Family::FIELDS)
        if (field.name != nullptr && other.name != nullptr &&
            std::string_view(field.name) == other.name)
          return true;
    return false;
  }
}

/**
 * @brief Whether @p Resource shares a field name with any resource of a list.
 * @tparam Resource Resource checked.
 * @tparam Resources Resource list.
 * @return True on any overlap.
 */
template <typename Resource, typename... Resources>
consteval bool resource_names_overlap(ResourceList<Resources...>) {
  return (resource_names_overlap<Resource, Resources>() || ... || false);
}

/**
 * @brief Whether a resource's slider names and noise seed carry its instance
 *        key.
 * @tparam Resource Parameter resource.
 * @tparam List The effect's resource list.
 * @return True for a non-standard key, a repeated family or a field-name
 *         overlap.
 */
template <typename Resource, typename List> consteval bool qualified() {
  return !Resource::STANDARD ||
         family_instances<typename Resource::Family>(List{}) > 1 ||
         resource_names_overlap<Resource>(List{});
}

/**
 * @brief Arena bytes the qualified slider names of one resource need.
 * @tparam Spec Composed-effect Spec, for gate evaluation.
 * @tparam Resource Parameter resource.
 * @tparam List The effect's resource list.
 * @return Bytes of "key.name" strings with terminators; 0 when unqualified or
 *         a colour family.
 */
template <typename Spec, typename Resource, typename List>
consteval size_t resource_name_bytes() {
  constexpr bool QUALIFY = qualified<Resource, List>();
  if constexpr (!QUALIFY || Resource::KIND == ResourceKind::COLOR)
    return 0;
  else {
    constexpr size_t PREFIX = Resource::KEY.view().size() + 2;
    if constexpr (Resource::KIND == ResourceKind::LENS) {
      size_t bytes = 0;
      for (const char *name : MOBIUS_PARAM_NAMES)
        bytes += PREFIX + std::string_view(name).size();
      return bytes;
    } else {
      size_t bytes = 0;
      for (const auto &field : Resource::Family::FIELDS)
        if (field.name != nullptr && field_gate_open<Spec>(field.gate))
          bytes += PREFIX + std::string_view(field.name).size();
      if constexpr (Resource::KIND == ResourceKind::WARP) {
        constexpr std::string_view SPEED_NAME =
            warp_speed_name(Resource::KEY.view());
        bytes += PREFIX + SPEED_NAME.size();
      }
      return bytes;
    }
  }
}
/**
 * @brief Arena bytes all qualified slider names of a resource list need.
 * @tparam Spec Composed-effect Spec.
 * @tparam Resources Parameter resources.
 * @return Sum of `resource_name_bytes` over @p Resources.
 */
template <typename Spec, typename... Resources>
consteval size_t parameter_name_bytes(ResourceList<Resources...>) {
  return (resource_name_bytes<Spec, Resources, ResourceList<Resources...>>() +
          ... + 0);
}

} // namespace ComposedDetail
