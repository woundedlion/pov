/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/** @file composed_descriptors.h
 * @brief Compile-time stage classification, topology derivation and parameter
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
                         typename Pipeline::template StageAt<Index>>::value)
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

/** @brief Semantic leaves without binding executable stages.
 * @tparam Pipeline Typed pipeline declaration. */
template <typename Pipeline> struct PipelineLeaves;
/** @brief Normalized leaf view of a typed pipeline.
 * @tparam B Declaration binding.
 * @tparam Entries Stage descriptors and placement groups. */
template <typename B, typename... Entries>
struct PipelineLeaves<Pullback::Pipeline<B, Entries...>> {
  using Nodes =
      typename Detail::NormalizeEntries<Entries...>::Type; ///< Execution nodes.
  using Leaves =
      typename Detail::Leaves<Nodes>::Type; ///< Flattened stage descriptors.
  static constexpr size_t STAGE_COUNT =
      Detail::LeafShape<Leaves>::COUNT; ///< Leaf count.
  /** @brief Leaf at an index. @tparam Index Leaf index. */
  template <size_t Index>
  using StageAt = typename Detail::TypeAt<Index, Leaves>::Type;
  /** @brief Whether any leaf satisfies a trait. @tparam Predicate Stage trait. */
  template <template <typename> class Predicate>
  static constexpr bool ANY_STAGE =
      Detail::LeafOps<Leaves>::template ANY<Predicate>;
  /** @brief First matching leaf, or void. @tparam Predicate Stage trait. */
  template <template <typename> class Predicate>
  using StageMatching =
      typename Detail::LeafOps<Leaves>::template Matching<Predicate>;
};

/** @brief Projection kind of a typed policy. @tparam Policy Projection policy. */
template <typename Policy> struct ProjectionKindOf {
  static_assert(!std::is_same_v<Policy, Policy>,
                "composed metadata: unsupported projection policy");
};
/** @brief Stereographic projection. @tparam P Projection provider. */
template <typename P> struct ProjectionKindOf<Projection::Stereographic<P>> {
  static constexpr auto VALUE =
      ProjectionKind::STEREOGRAPHIC; ///< Projection kind.
};
/** @brief Folded gnomonic projection. @tparam P Projection provider. */
template <typename P>
struct ProjectionKindOf<
    Projection::Gnomonic<P, Projection::GnomonicHemisphere::FOLDED>> {
  static constexpr auto VALUE =
      ProjectionKind::GNOMONIC_FOLDED; ///< Projection kind.
};
/** @brief Equirectangular projection. @tparam P Projection provider. */
template <typename P> struct ProjectionKindOf<Projection::Equirectangular<P>> {
  static constexpr auto VALUE =
      ProjectionKind::EQUIRECTANGULAR; ///< Projection kind.
};
/** @brief Folded sinusoidal projection. @tparam P Projection provider. */
template <typename P> struct ProjectionKindOf<Projection::FoldedSinusoidal<P>> {
  static constexpr auto VALUE =
      ProjectionKind::FOLDED_SINUSOIDAL; ///< Projection kind.
};
/** @brief Transfer kind of a typed policy. @tparam Policy Transfer policy. */
template <typename Policy> struct TransferKindOf {
  static_assert(!std::is_same_v<Policy, Policy>,
                "composed metadata: unsupported transfer policy");
};
/** @brief Iso-contour transfer. @tparam P Value provider. */
template <typename P> struct TransferKindOf<Transfer::IsoContour<P>> {
  static constexpr auto VALUE = TransferKind::ISO_CONTOUR; ///< Transfer kind.
};
/** @brief Field coverage kind of a typed policy. @tparam Policy Coverage policy. */
template <typename Policy> struct FieldCoverageKindOf {
  static_assert(!std::is_same_v<Policy, Policy>,
                "composed metadata: unsupported field coverage policy");
};
/** @brief Value cutout. @tparam P Value provider. */
template <typename P>
struct FieldCoverageKindOf<ValueCoverage::ValueCutout<P>> {
  static constexpr auto VALUE =
      FieldCoverageKind::VALUE_CUTOUT; ///< Coverage kind.
};
/** @brief Color modes encoded by a typed color policy. @tparam Policy Color policy. */
template <typename Policy> struct ColorModesOf {
  static_assert(!std::is_same_v<Policy, Policy>,
                "composed metadata: unsupported color policy");
};
/** @brief Generated-palette color modes.
 * @tparam B Declaration binding.
 * @tparam Hue Hue rotation mode.
 * @tparam Brightness Brightness envelope. */
template <typename B, HueMode Hue, Color::BrightnessEnvelope Brightness>
struct ColorModesOf<
    Color::GeneratedPalette<ColorProvider<B, Hue, Brightness>>> {
  static constexpr auto HUE = Hue;               ///< Hue rotation mode.
  static constexpr auto BRIGHTNESS = Brightness; ///< Brightness envelope.
};
/** @brief Whether a leaf colorizes sampled values. @tparam T Stage descriptor. */
template <typename T> struct IsColorStage : std::false_type {};
template <typename P>
struct IsColorStage<Stage::Colorize<P>> : std::true_type {};

/** @brief Topology derived from normalized typed pipeline entries.
 * @tparam Pipeline Semantic leaf view. */
template <typename Pipeline> struct PipelineMetadata {
  /// Whether a warp or surface tracks path length.
  static constexpr bool PATH_TRACKED =
      Pipeline::template ANY_STAGE<StagePathTracked>;
  using LensStage = typename Pipeline::template StageMatching<
      IsLensStage>; ///< First lens stage, or void.
  using ProjectStage = typename Pipeline::template StageMatching<
      IsProjectStage>; ///< Projection stage, or void.
  using SampleStage = typename Pipeline::template StageMatching<
      IsSampleStage>; ///< Sample stage, or void.
  using TransferStage = typename Pipeline::template StageMatching<
      IsTransferStage>; ///< Transfer stage, or void.
  using CoverageStage = typename Pipeline::template StageMatching<
      IsCoverageStage>; ///< Field coverage stage, or void.
  using ColorStage = typename Pipeline::template StageMatching<
      IsColorStage>; ///< Color stage, or void.
  static_assert(!std::is_void_v<ProjectStage>,
                "composed metadata: projection stage is required");
  static_assert(!std::is_void_v<ColorStage>,
                "composed metadata: color stage is required");
  /// Sphere-to-plane projection kind.
  static constexpr auto PROJECTION =
      ProjectionKindOf<typename ProjectStage::ProjectionPolicy>::VALUE;
  /// Projection coverage applied at sampling.
  static constexpr auto COVERAGE = [] {
    if constexpr (std::is_void_v<SampleStage>)
      return ProjectionCoverageMode::NONE;
    else
      return ProjectionCoverageModeOf<
          typename SampleStage::CoveragePolicy>::VALUE;
  }();
  static_assert(COVERAGE != static_cast<ProjectionCoverageMode>(255),
                "composed metadata: unsupported projection coverage policy");
  /// Sampled-value transfer kind.
  static constexpr auto TRANSFER = [] {
    if constexpr (std::is_void_v<TransferStage>)
      return TransferKind::NONE;
    else
      return TransferKindOf<typename TransferStage::TransferPolicy>::VALUE;
  }();
  /// Value-dependent field coverage kind.
  static constexpr auto FIELD_COVERAGE = [] {
    if constexpr (std::is_void_v<CoverageStage>)
      return FieldCoverageKind::NONE;
    else
      return FieldCoverageKindOf<typename CoverageStage::CoveragePolicy>::VALUE;
  }();
  /// Non-Mobius lens policy, or void.
  using LensPolicy = typename decltype([] {
    if constexpr (std::is_void_v<LensStage>)
      return std::type_identity<void>{};
    else if constexpr (IsMobiusLens<typename LensStage::LensPolicy>::value)
      return std::type_identity<void>{};
    else
      return std::type_identity<typename LensStage::LensPolicy>{};
  }())::type;
  /// Surface position relative to the lens; BEFORE_LENS when either is absent.
  static constexpr auto SURFACE_PLACEMENT = [] {
    constexpr size_t SURFACE = stage_index<Pipeline, IsSurfaceStage>();
    constexpr size_t LENS = stage_index<Pipeline, IsLensStage>();
    if constexpr (SURFACE == Pipeline::STAGE_COUNT ||
                  LENS == Pipeline::STAGE_COUNT || SURFACE < LENS)
      return SurfacePlacement::BEFORE_LENS;
    else
      return SurfacePlacement::AFTER_LENS;
  }();
  using ColorModes = ColorModesOf<
      typename ColorStage::ColorPolicy>;       ///< Typed color selections.
  static constexpr auto HUE = ColorModes::HUE; ///< Hue rotation mode.
  static constexpr auto BRIGHTNESS =
      ColorModes::BRIGHTNESS; ///< Brightness envelope.
};
/** @brief Topology derived without constructing a runtime frame.
 * @tparam Spec Authored pipeline spec. */
template <typename Spec>
struct MetadataFor
    : PipelineMetadata<
          PipelineLeaves<typename Spec::template Pipeline<DiscoveryBinding>>> {
  static_assert(
      !requires { Spec::PROJECTION; } && !requires { Spec::TRANSFER; } &&
          !requires { Spec::COVERAGE; } &&
          !requires { Spec::FIELD_COVERAGE; } && !requires { Spec::HUE; } &&
          !requires { Spec::BRIGHTNESS; } &&
          !requires { Spec::SURFACE_PLACEMENT; } &&
          !requires { typename Spec::LensPolicy; },
      "composed metadata: declare topology through typed pipeline entries");
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
 * @brief Whether a field gate is live under the derived projection topology.
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
    return uses_central_meridian(MetadataFor<Spec>::PROJECTION);
  case FieldGate::SINGULARITY_FADE:
    return uses_singularity_fade(MetadataFor<Spec>::PROJECTION);
  }
  return false;
}

/** Declared only: a call makes the enclosing constant evaluation ill-formed. */
void color_gate_names_unknown_topology();

/**
 * @brief Whether a colour field's gate is live under the typed color policy.
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
    value = static_cast<unsigned>(MetadataFor<Spec>::HUE);
  else if (topology == Color::BRIGHTNESS_ENVELOPE_GATE.field)
    value = static_cast<unsigned>(MetadataFor<Spec>::BRIGHTNESS);
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

/** @brief Topology inferred from a spec's typed pipeline. @tparam Spec Authored pipeline spec. */
template <typename Spec> using MetadataFor = ComposedDetail::MetadataFor<Spec>;
