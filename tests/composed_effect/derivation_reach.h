/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Derivation reach of the composed layer over the operator catalog
// ============================================================================

namespace In = Pullback::Interp;

using ReachBinding = Pullback::ComposedDetail::DiscoveryBinding;

template <typename Family>
concept DerivableSource = requires {
  typename Pullback::SourcePolicyFor<Family, ReachBinding>::Type;
};

template <typename Family>
concept DerivableWarp = requires {
  typename Pullback::WarpPolicyFor<Family, ReachBinding, "warp", false>::Type;
};

template <typename Family>
concept DerivableSurface = requires {
  typename Pullback::SurfacePolicyFor<Family, ReachBinding, false>::Type;
};

/** @brief The composed family an operator model's chain family extends.
 * @details The chain family adds catalog-only topology to it; the composed
 * layer fixes that topology in the derived policy.
 * @tparam Model Operator model. */
template <typename Model> struct ComposedFamily {
  using Type = typename Model::Params; ///< The family.
};
/** @brief Affine controls; the path-length flag is a policy argument. */
template <> struct ComposedFamily<In::Op::WarpAffineV4> {
  using Type = Pullback::AffineParams; ///< The family.
};
/** @brief Planar noise with the simplex basis. */
template <> struct ComposedFamily<In::Op::SampleProjectedNoise> {
  using Type = Pullback::ProjectedNoiseSourceParams; ///< The family.
};
/** @brief Ripple bursts; the chain family adds burst topology selectors. */
template <> struct ComposedFamily<In::Op::DisplaceRipple> {
  using Type = Pullback::SphereAnimation::BurstRippleParams; ///< The family.
};

/** @brief Lens operators whose composed policies a Spec names directly; no
 * family-keyed selector derives them. */
constexpr std::string_view SPEC_NAMED_LENSES[] = {
    "sphere.lens.glitch.v2",
    "sphere.lens.twist.v2",
    "sphere.lens.mobius.v3",
    "sphere.lens.kaleidoscope.v2",
};

/** @brief Whether every fixed tessellation kind and fill derives a source.
 * @tparam Kinds Kind values.
 * @return True when each kind derives under both fills. */
template <size_t... Kinds>
consteval bool tessellations_derivable(std::index_sequence<Kinds...>) {
  using Kind = Pattern::TessellationKind;
  using Fill = Pattern::TessellationFill;
  return ((DerivableSource<Pullback::FixedTessellationParams<
               static_cast<Kind>(Kinds), Fill::CELLS>> &&
           DerivableSource<Pullback::FixedTessellationParams<
               static_cast<Kind>(Kinds), Fill::EDGES>>) &&
          ...);
}

/** @brief Whether the composed layer derives a policy for an operator's
 * family: a source, warp or displacement, or a Spec-named lens.
 * @tparam Model Operator model.
 * @return True when derivable, or when the operator is not a source, warp or
 *         lens. */
template <typename Model> consteval bool family_derivable() {
  using Family = typename ComposedFamily<Model>::Type;
  constexpr In::CarrierId INPUT = In::carrier_id_of<typename Model::Input>();
  constexpr In::CarrierId OUTPUT = In::carrier_id_of<typename Model::Output>();
  if constexpr (std::is_same_v<Model, In::Op::SampleTessellation>)
    return tessellations_derivable(
        std::make_index_sequence<std::size(In::Op::TESSELLATION_KIND_IDS)>{});
  else if constexpr (OUTPUT == In::CarrierId::FIELD &&
                     INPUT != In::CarrierId::FIELD)
    return DerivableSource<Family>;
  else if constexpr (INPUT == In::CarrierId::PLANE &&
                     OUTPUT == In::CarrierId::PLANE)
    return DerivableWarp<Family>;
  else if constexpr (INPUT == In::CarrierId::SPHERE &&
                     OUTPUT == In::CarrierId::SPHERE &&
                     !std::is_same_v<Model, In::Op::Rotate>)
    return DerivableSurface<Family> ||
           std::ranges::find(SPEC_NAMED_LENSES, std::string_view(Model::ID)) !=
               std::end(SPEC_NAMED_LENSES);
  else
    return true;
}

/** @brief Whether every operator family derives a composed policy.
 * @tparam Models Operator models.
 * @return True when each model's family derives. */
template <typename... Models>
consteval bool families_derivable(Pullback::Detail::TypeList<Models...>) {
  return (family_derivable<Models>() && ...);
}

static_assert(families_derivable(In::OperatorModels{}));

// The shared noise field group is not itself a family, so it names no single
// policy.
static_assert(!DerivableSource<Pullback::Source::NoiseSourceParams>);

static_assert(
    std::is_same_v<
        typename Pullback::SourcePolicyFor<Pullback::ProjectedNoiseSourceParams,
                                           ReachBinding>::Type,
        Pullback::Source::ProjectedNoise<
            Pullback::SourceProvider<ReachBinding,
                                     Pullback::ProjectedNoiseSourceParams>,
            math::NoiseBasis::SIMPLEX>>);
static_assert(std::is_same_v<
              typename Pullback::SourcePolicyFor<
                  Pullback::FixedTessellationParams<
                      Pattern::TessellationKind::HEXAGONAL,
                      Pattern::TessellationFill::EDGES>,
                  ReachBinding>::Type,
              Pullback::Source::Tessellation<
                  Pullback::SourceProvider<
                      ReachBinding, Pullback::FixedTessellationParams<
                                        Pattern::TessellationKind::HEXAGONAL,
                                        Pattern::TessellationFill::EDGES>>,
                  Pattern::TessellationKind::HEXAGONAL,
                  Pattern::TessellationFill::EDGES>>);
static_assert(
    std::is_same_v<
        typename Pullback::WarpPolicyFor<Pullback::PolarParams, ReachBinding,
                                         "outer_warp", false>::Type,
        Pullback::Warp::PolarChart<
            Pullback::WarpProvider<ReachBinding, "outer_warp",
                                   Pullback::PolarParams, false>,
            Pullback::Warp::LinearPolar, 1>>);
static_assert(std::is_same_v<
              typename Pullback::WarpPolicyFor<Pullback::VectorNoiseParams,
                                               ReachBinding, "inner_warp",
                                               false>::Type,
              Pullback::Warp::VectorNoise<
                  Pullback::WarpProvider<ReachBinding, "inner_warp",
                                         Pullback::VectorNoiseParams, false>,
                  math::NoiseBasis::SIMPLEX, Pullback::Warp::FlatEnvelope>>);
static_assert(
    std::is_same_v<typename Pullback::SurfacePolicyFor<
                       Pullback::SurfaceNoiseParams, ReachBinding, false>::Type,
                   Pullback::Surface::CurlNoise<
                       Pullback::SurfaceProvider<
                           ReachBinding, Pullback::SurfaceNoiseParams, false>,
                       math::NoiseBasis::SIMPLEX, Pullback::Surface::Euler,
                       math::TangentLimit::SMOOTH>>);
static_assert(
    std::is_same_v<
        typename Pullback::SurfacePolicyFor<
            Pullback::AuthoredSurfaceNoiseParams, ReachBinding, false>::Type,
        Pullback::Surface::ConfiguredCurl<Pullback::SurfaceProvider<
            ReachBinding, Pullback::AuthoredSurfaceNoiseParams, false>>>);
