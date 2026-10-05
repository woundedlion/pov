/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

namespace ComposedDetail {

inline constexpr const char *MOBIUS_PARAM_NAMES[] = {
    "Mobius A Re", "Mobius A Im", "Mobius B Re", "Mobius B Im",
    "Mobius C Re", "Mobius C Im", "Mobius D Re", "Mobius D Im"};
constexpr std::string_view warp_speed_name(std::string_view key) {
  return key == "outer_warp"   ? "Planar Warp 1 Speed"
         : key == "inner_warp" ? "Planar Warp 2 Speed"
                               : "Planar Warp Speed";
}

template <typename T> struct IsSampleStage : std::false_type {};
template <typename S, typename W, typename C>
struct IsSampleStage<Stage::Sample<S, W, C>> : std::true_type {};
template <typename T> struct ProjectionCoverageModeOf {
  static constexpr auto VALUE = static_cast<ProjectionCoverageMode>(255);
};
template <> struct ProjectionCoverageModeOf<ProjectionCoverage::None> {
  static constexpr auto VALUE = ProjectionCoverageMode::NONE;
};
template <> struct ProjectionCoverageModeOf<ProjectionCoverage::Weight> {
  static constexpr auto VALUE = ProjectionCoverageMode::WEIGHT;
};
template <> struct ProjectionCoverageModeOf<ProjectionCoverage::WeightSquared> {
  static constexpr auto VALUE = ProjectionCoverageMode::WEIGHT_SQUARED;
};
template <typename P>
struct ProjectionCoverageModeOf<ProjectionCoverage::EdgeFade<P>> {
  static constexpr auto VALUE = ProjectionCoverageMode::EDGE_FADE;
};

template <typename T> struct IsLensStage : std::false_type {};
template <typename P> struct IsLensStage<Stage::Lens<P>> : std::true_type {};
template <typename T> struct IsSurfaceStage : std::false_type {};
template <typename P>
struct IsSurfaceStage<Stage::Displace<P>> : std::true_type {};
template <typename T> struct IsProjectStage : std::false_type {};
template <typename P>
struct IsProjectStage<Stage::Project<P>> : std::true_type {};
template <typename T> struct IsTransferStage : std::false_type {};
template <typename P>
struct IsTransferStage<Stage::Transfer<P>> : std::true_type {};
template <typename T> struct IsCoverageStage : std::false_type {};
template <typename P>
struct IsCoverageStage<Stage::ApplyCoverage<P>> : std::true_type {};
template <typename T> struct IsMobiusLens : std::false_type {};
template <typename P> struct IsMobiusLens<Lens::Mobius<P>> : std::true_type {};

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
template <typename Provider, math::NoiseBasis Basis, typename Integrator>
struct PathTracked<Surface::CurlNoise<Provider, Basis, Integrator>>
    : PathTracked<Provider> {};
template <typename Provider, math::NoiseBasis Basis>
struct PathTracked<Surface::DirectNoise<Provider, Basis>>
    : PathTracked<Provider> {};
template <typename Stage>
struct StagePathTracked : PathTracked<typename Stage::Policies> {};

template <typename Spec, typename Binding, typename Pipeline>
struct PipelineMetadata {
  static constexpr bool PATH_TRACKED =
      Pipeline::template any_stage<StagePathTracked>;
  using LensStage = typename Pipeline::template stage_matching<IsLensStage>;
  using ProjectStage =
      typename Pipeline::template stage_matching<IsProjectStage>;
  using SampleStage = typename Pipeline::template stage_matching<IsSampleStage>;
  static constexpr bool COVERAGE_MATCHES = [] {
    if constexpr (std::is_void_v<SampleStage>)
      return Spec::COVERAGE == ProjectionCoverageMode::NONE;
    else
      return ProjectionCoverageModeOf<
                 typename SampleStage::CoveragePolicy>::VALUE == Spec::COVERAGE;
  }();
  static constexpr bool LENS_MATCHES = [] {
    if constexpr (std::is_void_v<LensStage>)
      return std::is_void_v<typename Spec::LensPolicy>;
    else if constexpr (IsMobiusLens<typename LensStage::LensPolicy>::value)
      return std::is_void_v<typename Spec::LensPolicy>;
    else
      return std::is_same_v<typename LensStage::LensPolicy,
                            typename Spec::LensPolicy>;
  }();
  static constexpr bool PROJECTION_MATCHES = [] {
    if constexpr (std::is_void_v<ProjectStage>)
      return false;
    else
      return std::is_same_v<
          typename ProjectStage::ProjectionPolicy,
          typename ProjectionPolicyFor<Spec::PROJECTION, Binding>::Type>;
  }();
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

template <typename B> struct PolicyResources<OuterCameraProvider<B>> {
  using Type = ResourceList<ParameterResource<"projection", ProjectionParams,
                                              ResourceKind::PROJECTION>>;
};
template <typename B>
struct PolicyResources<ProjectionProvider<B>>
    : PolicyResources<OuterCameraProvider<B>> {};
template <typename B, ResourceKey Key>
struct PolicyResources<LensProvider<B, Key>> {
  using Type = ResourceList<
      ParameterResource<Key, MobiusLensParams, ResourceKind::LENS>>;
};
template <typename B, ResourceKey Key, typename Family, bool Track,
          ResourceKey SourceKey>
struct PolicyResources<WarpProvider<B, Key, Family, Track, SourceKey>> {
  using Type = ResourceList<ParameterResource<Key, Family, ResourceKind::WARP>>;
};
template <typename B, typename Family, bool Track, ResourceKey Key>
struct PolicyResources<SurfaceProvider<B, Family, Track, Key>> {
  using Type =
      ResourceList<ParameterResource<Key, Family, ResourceKind::SURFACE>>;
};
template <typename B, typename Family, ResourceKey Key>
struct PolicyResources<SourceProvider<B, Family, Key>> {
  using Type =
      ResourceList<ParameterResource<Key, Family, ResourceKind::SOURCE>>;
};
template <typename B, typename Family, ResourceKey Key>
struct PolicyResources<ValueProvider<B, Family, Key>> {
  using Type =
      ResourceList<ParameterResource<Key, Family, ResourceKind::VALUE>>;
};
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
template <typename Provider, math::NoiseBasis Basis, typename Integrator>
struct PolicyResources<Surface::CurlNoise<Provider, Basis, Integrator>>
    : PolicyResources<Provider> {};
template <typename Provider, math::NoiseBasis Basis>
struct PolicyResources<Surface::DirectNoise<Provider, Basis>>
    : PolicyResources<Provider> {};

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
        if (field.member == &ColorParams::hue_shift_amount &&
            Spec::HUE == HueMode::NONE)
          continue;
        if ((field.member == &ColorParams::hue_noise_scale ||
             field.member == &ColorParams::hue_noise_speed) &&
            Spec::HUE != HueMode::NOISE)
          continue;
        if ((field.member == &ColorParams::brightness_bottom ||
             field.member == &ColorParams::brightness_top) &&
            Spec::BRIGHTNESS == Color::BrightnessEnvelope::NONE)
          continue;
        ++count;
      } else if (field.name != nullptr && field_gate_open<Spec>(field.gate))
        ++count;
    }
    return count;
  }
}
template <typename Spec, typename... Resources>
consteval size_t parameter_count(ResourceList<Resources...>) {
  return (resource_parameter_count<Spec, Resources>() + ... + 0);
}

template <typename Family, typename... Resources>
consteval size_t family_instances(ResourceList<Resources...>) {
  return (size_t(std::is_same_v<Family, typename Resources::Family>) + ... + 0);
}

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

template <typename Resource, typename... Resources>
consteval bool resource_names_overlap(ResourceList<Resources...>) {
  return (resource_names_overlap<Resource, Resources>() || ... || false);
}

template <typename Resource, typename List> consteval bool qualified() {
  return !Resource::STANDARD ||
         family_instances<typename Resource::Family>(List{}) > 1 ||
         resource_names_overlap<Resource>(List{});
}

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
template <typename Spec, typename... Resources>
consteval size_t parameter_name_bytes(ResourceList<Resources...>) {
  return (resource_name_bytes<Spec, Resources, ResourceList<Resources...>>() +
          ... + 0);
}

} // namespace ComposedDetail
