/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/** @file composed_policies.h
 * @brief Composed-effect authored controls and policy selectors keyed on
 * parameter families and projection kind.
 */

/** @brief Authored controls accompanying an explicit typed stage pipeline.
 * @details Derived specs supply `template <typename B> using Pipeline`.
 * Stage topology is exposed by `MetadataFor`. Proven performance paths are on
 * by default; an effect sets one false only for a constraint it names. */
struct Spec {
  static constexpr bool AUTHORED_TOPOLOGY =
      false; ///< Uses canonical authored node identities.
  /// Hue harmony of generated palettes.
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  /// Whether the projection frame random-walks and spins.
  static constexpr bool ANIMATED_PROJECTION = true;
  /// Whether columns sample consecutive frame states across a half-turn.
  static constexpr bool TEMPORAL_SAMPLING = false;
  /// Samples eligible color-source phases once after palette interpolation;
  /// effective only with `TEMPORAL_SAMPLING` and without fusion.
  static constexpr bool TEMPORAL_COLOR_PHASE = true;
  /// Fuses an eligible final temporal color layer.
  static constexpr bool FUSED_TEMPORAL_COLOR = false;
  /// Bakes palette-coordinate by hue-rotation tables.
  static constexpr bool HUE_ROTATION_LUT = true;
  /// Caches RGB-first interval palette blends in OKLab; effective only with
  /// column phase sampling or fused temporal color.
  static constexpr bool PALETTE_LAB_CACHE = true;
  /// Caches palette colors after hue rotation and gamut mapping; effective
  /// only with a palette LAB cache on the hue-rotation table's value grid and
  /// the NONE brightness envelope.
  static constexpr bool PALETTE_RGB_CACHE = true;
  /// Caches direct noise over the dodecahedral chamber at temporal endpoints.
  static constexpr bool PHASE_NOISE_CACHE = false;
  static constexpr bool PALETTE_RGB_MIDPOINT_ONLY =
      false; ///< Stores one admitted RGB midpoint.
  /// Limits interval caches to the initial color and modulation controls.
  static constexpr bool INITIAL_CACHE_CONTROLS_ONLY = false;
  /// Palette-coordinate nodes per cached time row.
  static constexpr int PALETTE_LAB_VALUES = 64;
  /// Factors a supported rotation prefix into column bases.
  static constexpr bool FACTORED_COLUMNS = false;
  /// Bakes endpoint radial transfer curves for burst ripples.
  static constexpr bool RIPPLE_RADIAL_TABLE = false;
  /// Keeps each prepared shading frame in the runtime arena.
  static constexpr bool PERSISTENT_PREPARED_FRAME = true;
  /// Keeps the full-resolution gamut table in device persistent storage.
  static constexpr bool CACHE_FULL_GAMUT_LUT = true;
  /// Copies immutable dodecahedral reflection matrices into persistent storage.
  static constexpr bool CACHE_DODECAHEDRAL_MATRICES = false;
  static constexpr bool CACHE_DODECAHEDRAL_FACES =
      false; ///< Copies cube classifiers into persistent storage.
};

/** @brief Catalog identity of one projection algorithm.
 * @param projection Projection algorithm.
 * @return Versioned operator identity, or empty for no projection. */
constexpr std::string_view projection_operator_id(ProjectionKind projection) {
  switch (projection) {
  case ProjectionKind::STEREOGRAPHIC:
    return "project.stereographic.v2";
  case ProjectionKind::GNOMONIC_FOLDED:
  case ProjectionKind::GNOMONIC:
    return "project.gnomonic.v2";
  case ProjectionKind::EQUIRECTANGULAR:
    return "project.equirectangular.v2";
  case ProjectionKind::FOLDED_SINUSOIDAL:
    return "project.folded-sinusoidal.v2";
  case ProjectionKind::PEIRCE_SQUARE_FAST:
    return "project.peirce-square-fast.v3";
  case ProjectionKind::PEIRCE:
    return "project.peirce.v3";
  case ProjectionKind::BONNE:
    return "project.bonne.v3";
  case ProjectionKind::AIROCEAN:
    return "project.airocean.v3";
  case ProjectionKind::TETRAHEDRAL:
    return "project.tetrahedral.v1";
  case ProjectionKind::CUBIC:
    return "project.cubic.v1";
  case ProjectionKind::OCTAHEDRAL:
    return "project.octahedral.v1";
  case ProjectionKind::DODECAHEDRAL:
    return "project.dodecahedral.v1";
  case ProjectionKind::ICOSAHEDRAL:
    return "project.icosahedral.v1";
  case ProjectionKind::NONE:
    return {};
  }
  return {};
}

/**
 * @brief Maps a source parameter family to its `Pullback::Source` policy,
 *        exposed as `Type`.
 * @tparam Family Source parameter family.
 * @tparam Binding Frame binding of the effect.
 */
template <typename Family, typename Binding> struct SourcePolicyFor;
/** @brief Grid source. @tparam B Frame binding. */
template <typename B> struct SourcePolicyFor<GridSourceParams, B> {
  using Type = Pullback::Source::Grid<SourceProvider<B, GridSourceParams>>;
};
/** @brief Twin-wave source. @tparam B Frame binding. */
template <typename B> struct SourcePolicyFor<TwinWaveSourceParams, B> {
  using Type =
      Pullback::Source::TwinWave<SourceProvider<B, TwinWaveSourceParams>>;
};
/** @brief Plane bands source. @tparam B Frame binding. */
template <typename B> struct SourcePolicyFor<BandsSourceParams, B> {
  using Type = Pullback::Source::Bands<SourceProvider<B, BandsSourceParams>>;
};
/** @brief Spiral source. @tparam B Frame binding. */
template <typename B> struct SourcePolicyFor<SpiralSourceParams, B> {
  using Type = Pullback::Source::Spiral<SourceProvider<B, SpiralSourceParams>>;
};
/** @brief Fixed-kind tessellation source.
 * @tparam Kind Cell shape.
 * @tparam Fill Cell primitives or edge band.
 * @tparam B Frame binding. */
template <Pattern::TessellationKind Kind, Pattern::TessellationFill Fill,
          typename B>
struct SourcePolicyFor<FixedTessellationParams<Kind, Fill>, B> {
  using Type = Pullback::Source::Tessellation<
      SourceProvider<B, FixedTessellationParams<Kind, Fill>>, Kind, Fill>;
};
/** @brief Planar simplex-noise source. @tparam B Frame binding. */
template <typename B> struct SourcePolicyFor<ProjectedNoiseSourceParams, B> {
  using Type = Pullback::Source::ProjectedNoise<
      SourceProvider<B, ProjectedNoiseSourceParams>, math::NoiseBasis::SIMPLEX>;
};

/** @brief Authored source topology. @tparam B Frame binding. */
template <typename B>
struct SourcePolicyFor<Source::AuthoredProjectedNoiseSourceParams, B> {
  using Type = Source::ConfiguredProjectedNoise<
      SourceProvider<B, Source::AuthoredProjectedNoiseSourceParams>>;
};

namespace ComposedDetail {
/** @brief Rebinds a source policy to authored crossing controls.
 * @tparam Policy Existing concrete source policy.
 * @tparam Provider Authored crossing provider. */
template <typename Policy, typename Provider> struct RebindSource;
/** @brief Single-provider source policy.
 * @tparam Policy Source policy template.
 * @tparam Old Original provider.
 * @tparam Provider Authored crossing provider. */
template <template <typename> typename Policy, typename Old, typename Provider>
struct RebindSource<Policy<Old>, Provider> {
  using Type = Policy<Provider>;
};
/** @brief Fixed-kind tessellation policy.
 * @tparam Old Original provider.
 * @tparam Kind Cell shape.
 * @tparam Fill Fill.
 * @tparam Provider Authored crossing provider. */
template <typename Old, Pattern::TessellationKind Kind,
          Pattern::TessellationFill Fill, typename Provider>
struct RebindSource<Source::Tessellation<Old, Kind, Fill>, Provider> {
  using Type = Source::Tessellation<Provider, Kind, Fill>;
};
} // namespace ComposedDetail

/** @brief Twin-wave source with the full authored drift range. @tparam B Frame binding. */
template <typename B>
struct SourcePolicyFor<Source::AuthoredTwinWaveSourceParams, B> {
  using Type =
      Source::TwinWave<SourceProvider<B, Source::AuthoredTwinWaveSourceParams>>;
};
/** @brief Source policy retaining runtime crossing controls.
 * @tparam Family Source-specific family.
 * @tparam B Frame binding. */
template <typename Family, typename B>
struct SourcePolicyFor<Source::AuthoredSampleParams<Family>, B> {
  using Type = typename ComposedDetail::RebindSource<
      typename SourcePolicyFor<Family, B>::Type,
      SourceProvider<B, Source::AuthoredSampleParams<Family>>>::Type;
};

/**
 * @brief Maps a warp parameter family to its `Pullback::Warp` policy, exposed
 *        as `Type`.
 * @tparam Family Warp parameter family.
 * @tparam Binding Frame binding of the effect.
 * @tparam Key Parameter resource key of this warp instance.
 * @tparam TrackPath Whether the warp accumulates path length for hue.
 * @tparam Form Parts an affine warp evaluates; other families take the
 *         default.
 */
template <typename Family, typename Binding, ResourceKey Key, bool TrackPath,
          Warp::AffineForm Form = Warp::AffineForm::GENERAL>
struct WarpPolicyFor;
/**
 * @brief Mirror-tile warp.
 * @tparam B Frame binding.
 * @tparam K Parameter resource key.
 * @tparam T Whether path length is tracked.
 */
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<MirrorParams, B, K, T> {
  using Type = Pullback::Warp::MirrorTile<WarpProvider<B, K, MirrorParams, T>>;
};
/**
 * @brief Wave-shear warp.
 * @tparam B Frame binding.
 * @tparam K Parameter resource key.
 * @tparam T Whether path length is tracked.
 */
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<WaveShearParams, B, K, T> {
  using Type =
      Pullback::Warp::WaveShear<WarpProvider<B, K, WaveShearParams, T>>;
};
/**
 * @brief Simplex vector-noise warp, flat envelope.
 * @tparam B Frame binding.
 * @tparam K Parameter resource key.
 * @tparam T Whether path length is tracked.
 */
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<VectorNoiseParams, B, K, T> {
  using Type =
      Pullback::Warp::VectorNoise<WarpProvider<B, K, VectorNoiseParams, T>,
                                  math::NoiseBasis::SIMPLEX,
                                  Pullback::Warp::FlatEnvelope>;
};
/**
 * @brief Affine-frame warp.
 * @tparam B Frame binding.
 * @tparam K Parameter resource key.
 * @tparam T Whether path length is tracked.
 * @tparam F Parts the warp evaluates.
 */
template <typename B, ResourceKey K, bool T, Warp::AffineForm F>
struct WarpPolicyFor<AffineParams, B, K, T, F> {
  using Type = Pullback::Warp::AffineFrame<
      WarpProvider<B, K, AffineParams, T, "source", F>>;
};
/** @brief Orbiting vortex warp.
 * @tparam B Frame binding.
 * @tparam K Resource key.
 * @tparam T Whether path length is tracked. */
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<Warp::VortexParams, B, K, T> {
  using Type =
      Pullback::Warp::Vortex<WarpProvider<B, K, Warp::VortexParams, T>>;
};
/** @brief Euler simplex curl flow.
 * @tparam B Frame binding.
 * @tparam K Resource key.
 * @tparam T Whether path length is tracked. */
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<CurlFlowParams, B, K, T> {
  using Type = Pullback::Warp::CurlFlow<WarpProvider<B, K, CurlFlowParams, T>,
                                        math::NoiseBasis::SIMPLEX,
                                        Pullback::Warp::Euler1>;
};

/**
 * @brief First-harmonic linear polar-chart warp.
 * @tparam B Frame binding.
 * @tparam K Parameter resource key.
 * @tparam T Whether path length is tracked.
 */
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<PolarParams, B, K, T> {
  using Type = Pullback::Warp::PolarChart<WarpProvider<B, K, PolarParams, T>,
                                          Pullback::Warp::LinearPolar, 1>;
};

/** @brief Authored warp topology.
 * @tparam B Frame binding.
 * @tparam K Resource key.
 * @tparam T Whether path length is tracked. */
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<Warp::AuthoredWaveShearParams, B, K, T> {
  using Type = Warp::ConfiguredWaveShear<
      WarpProvider<B, K, Warp::AuthoredWaveShearParams, T>>;
};
/** @brief Authored warp topology.
 * @tparam B Frame binding.
 * @tparam K Resource key.
 * @tparam T Whether path length is tracked. */
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<Warp::AuthoredVectorNoiseParams, B, K, T> {
  using Type = Warp::ConfiguredVectorNoise<
      WarpProvider<B, K, Warp::AuthoredVectorNoiseParams, T>>;
};
/** @brief Authored warp topology.
 * @tparam B Frame binding.
 * @tparam K Resource key.
 * @tparam T Whether path length is tracked. */
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<Warp::AuthoredPolarParams, B, K, T> {
  using Type = Warp::ConfiguredPolarChart<
      WarpProvider<B, K, Warp::AuthoredPolarParams, T>>;
};
/** @brief Authored warp topology.
 * @tparam B Frame binding.
 * @tparam K Resource key.
 * @tparam T Whether path length is tracked. */
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<Warp::AuthoredCurlFlowParams, B, K, T> {
  using Type = Warp::ConfiguredCurlFlow<
      WarpProvider<B, K, Warp::AuthoredCurlFlowParams, T>>;
};

/**
 * @brief Maps a surface parameter family to its `Pullback::Surface` policy,
 *        exposed as `Type`.
 * @tparam Family Surface parameter family.
 * @tparam Binding Frame binding of the effect.
 * @tparam TrackPath Whether the surface accumulates path length for hue.
 */
template <typename Family, typename Binding, bool TrackPath>
struct SurfacePolicyFor;
/**
 * @brief Euler-integrated simplex curl-noise flow.
 * @tparam B Frame binding.
 * @tparam T Whether path length is tracked.
 */
template <typename B, bool T>
struct SurfacePolicyFor<SurfaceNoiseParams, B, T> {
  using Type = Pullback::Surface::CurlNoise<
      SurfaceProvider<B, SurfaceNoiseParams, T>, math::NoiseBasis::SIMPLEX,
      Pullback::Surface::Euler, math::TangentLimit::SMOOTH>;
};
/**
 * @brief Direct simplex-noise displacement.
 * @tparam B Frame binding.
 * @tparam T Whether path length is tracked.
 */
template <typename B, bool T>
struct SurfacePolicyFor<DirectSurfaceParams, B, T> {
  using Type =
      Pullback::Surface::DirectNoise<SurfaceProvider<B, DirectSurfaceParams, T>,
                                     math::NoiseBasis::SIMPLEX>;
};
/**
 * @brief Periodic ripple surface.
 * @tparam B Frame binding.
 * @tparam T Whether path length is tracked.
 */
template <typename B, bool T>
struct SurfacePolicyFor<PeriodicRippleParams, B, T> {
  using Type = Pullback::Surface::PeriodicRipple<
      SurfaceProvider<B, PeriodicRippleParams, T>>;
};

/**
 * @brief Maps a `ProjectionKind` to its `Pullback::Projection` policy, exposed
 *        as `Type`.
 * @tparam ProjectionV Projection kind.
 * @tparam Binding Frame binding of the effect.
 */
template <ProjectionKind ProjectionV, typename Binding>
struct ProjectionPolicyFor;
/** @brief Stereographic projection. @tparam B Frame binding. */
template <typename B>
struct ProjectionPolicyFor<ProjectionKind::STEREOGRAPHIC, B> {
  using Type = Pullback::Projection::Stereographic<ProjectionProvider<B>>;
};
/** @brief Gnomonic projection, hemispheres folded. @tparam B Frame binding. */
template <typename B>
struct ProjectionPolicyFor<ProjectionKind::GNOMONIC_FOLDED, B> {
  using Type = Pullback::Projection::Gnomonic<
      ProjectionProvider<B>, Pullback::Projection::GnomonicHemisphere::FOLDED>;
};
/** @brief Equirectangular projection. @tparam B Frame binding. */
template <typename B>
struct ProjectionPolicyFor<ProjectionKind::EQUIRECTANGULAR, B> {
  using Type = Pullback::Projection::Equirectangular<ProjectionProvider<B>>;
};
/** @brief Folded sinusoidal projection. @tparam B Frame binding. */
template <typename B>
struct ProjectionPolicyFor<ProjectionKind::FOLDED_SINUSOIDAL, B> {
  using Type = Pullback::Projection::FoldedSinusoidal<ProjectionProvider<B>>;
};

/** @brief Solid table of a polyhedral atlas projection, exposed as `Type`.
 * @tparam Kind Polyhedral atlas projection. */
template <ProjectionKind Kind> struct AtlasSolidFor;
/** @brief Tetrahedral atlas solid. */
template <> struct AtlasSolidFor<ProjectionKind::TETRAHEDRAL> {
  using Type = Solids::Tetrahedron; ///< Solid table.
};
/** @brief Cubic atlas solid. */
template <> struct AtlasSolidFor<ProjectionKind::CUBIC> {
  using Type = Solids::Cube; ///< Solid table.
};
/** @brief Octahedral atlas solid. */
template <> struct AtlasSolidFor<ProjectionKind::OCTAHEDRAL> {
  using Type = Solids::Octahedron; ///< Solid table.
};
/** @brief Dodecahedral atlas solid. */
template <> struct AtlasSolidFor<ProjectionKind::DODECAHEDRAL> {
  using Type = Solids::Dodecahedron; ///< Solid table.
};
/** @brief Icosahedral atlas solid. */
template <> struct AtlasSolidFor<ProjectionKind::ICOSAHEDRAL> {
  using Type = Solids::Icosahedron; ///< Solid table.
};

/** @brief Runtime topology for authored atlas projections.
 * @tparam Kind Projection algorithm.
 * @tparam State Parameter provider. */
template <ProjectionKind Kind, typename State>
struct AuthoredProjection
    : Projection::FrameConjugateFromState<State>,
      Projection::MeridianPrepared<AuthoredProjection<Kind, State>, State> {
  static_assert(Kind < ProjectionKind::NONE,
                "authored projection requires a supported projection kind");
  using Binding = typename State::Binding;       ///< Pipeline binding.
  using FrameState = typename State::FrameState; ///< Bound frame.
  using Prepared = Projection::PreparedMeridian; ///< Cached meridian rotation.
  using Projection::MeridianPrepared<AuthoredProjection, State>::project;
  static constexpr bool EDGE_DISTANCE_AVAILABLE =
      Kind != ProjectionKind::FOLDED_SINUSOIDAL; ///< Cut distances available.
  /// Whether the projection uses the approximate Peirce kernel.
  static constexpr bool APPROXIMATE =
      Kind == ProjectionKind::PEIRCE_SQUARE_FAST;
  /// Reference implementation for the approximate kernel.
  static constexpr ApproximationOracleId ORACLE =
      APPROXIMATE ? ApproximationOracleId::PEIRCE_FAST_SQUARE
                  : ApproximationOracleId::NONE;
  static constexpr auto METRICS = [] {
    if constexpr (APPROXIMATE)
      return Projection::PEIRCE_FAST_SQUARE_METRICS;
    else
      return std::array<ApproximationMetric, 0>{};
  }(); ///< Approximation error bounds.
  /** @brief Selects the authored projection-frame behavior.
   * @param frame Bound frame.
   * @return Identity or animated conjugate rotation. */
  static const math::Quaternion &frame_conjugate(const FrameState &frame) {
    static const math::Quaternion IDENTITY;
    return frame.params.template get<"projection">().frame == 0
               ? IDENTITY
               : State::conjugate(frame);
  }
  /** @brief Checks the bound projection provider.
   * @tparam CandidateBinding Candidate pipeline binding. */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Projection::FrameProvider<State, CandidateBinding>;
  /** @brief Projects with authored layout and hemisphere controls.
   * @param input Projection-local unit direction.
   * @param frame Bound frame.
   * @param prepared Cached meridian rotation.
   * @return Plane coordinates and projection provenance. */
  static ProjectionResult project(const math::Vector &input,
                                  const FrameState &frame,
                                  const Prepared &prepared) {
    const auto &params = frame.params.template get<"projection">();
    if constexpr (Kind == ProjectionKind::STEREOGRAPHIC)
      return Projection::stereographic(input, params.singularity_fade);
    else if constexpr (Kind == ProjectionKind::FOLDED_SINUSOIDAL)
      return Projection::folded_sinusoidal(input, params.central_meridian);
    else if constexpr (Kind == ProjectionKind::EQUIRECTANGULAR)
      return Projection::equirectangular(input, params.central_meridian,
                                         params.singularity_fade);
    else if constexpr (Kind == ProjectionKind::GNOMONIC_FOLDED)
      return Projection::gnomonic(
          input, params.singularity_fade,
          static_cast<Projection::GnomonicHemisphere>(params.hemisphere));
    else if constexpr (Kind == ProjectionKind::PEIRCE)
      return Projection::peirce(
          input, params.central_meridian, params.layout, params.layout_scroll,
          true, params.coordinate_scale, params.singularity_fade,
          prepared.cosine, prepared.sine);
    else if constexpr (Kind == ProjectionKind::PEIRCE_SQUARE_FAST)
      return Projection::peirce_fast_square(input, params.coordinate_scale,
                                            params.singularity_fade);
    else if constexpr (Kind == ProjectionKind::BONNE)
      return Projection::bonne(input, params.central_meridian,
                               (params.hemisphere == 0 ? 1.0f : -1.0f) *
                                   params.standard_parallel,
                               params.coordinate_scale);
    else if constexpr (Kind == ProjectionKind::AIROCEAN)
      return Projection::airocean(input, params.layout == 1, true,
                                  params.coordinate_scale, prepared.cosine,
                                  prepared.sine);
    else if constexpr (is_polyhedral_atlas(Kind))
      return Projection::from_kernel(
          projections::polyhedral_atlas_projection(
              projections::POLYHEDRAL_ATLAS<typename AtlasSolidFor<Kind>::Type>,
              input, prepared.cosine, prepared.sine, params.area_correction),
          params.coordinate_scale);
    else
      return Projection::gnomonic(
          input, params.singularity_fade,
          static_cast<Projection::GnomonicHemisphere>(params.hemisphere));
  }
};

/** @brief Selectable authored projection policy.
 * @tparam ProjectionV Projection algorithm.
 * @tparam B Pipeline binding. */
template <ProjectionKind ProjectionV, typename B>
  requires(ProjectionV == ProjectionKind::PEIRCE ||
           ProjectionV == ProjectionKind::PEIRCE_SQUARE_FAST ||
           ProjectionV == ProjectionKind::BONNE ||
           ProjectionV == ProjectionKind::AIROCEAN ||
           is_polyhedral_atlas(ProjectionV) ||
           ProjectionV == ProjectionKind::GNOMONIC)
struct ProjectionPolicyFor<ProjectionV, B> {
  using Type = AuthoredProjection<ProjectionV,
                                  ProjectionProvider<B, AtlasProjectionParams>>;
};

/**
 * @brief Maps a `ProjectionCoverageMode` to its `Pullback::ProjectionCoverage`
 *        policy, exposed as `Type`.
 * @tparam CoverageV Coverage mode.
 * @tparam Binding Frame binding of the effect.
 * @tparam Family Value parameter family; void when none is read.
 * @tparam Key Value resource key.
 */
template <ProjectionCoverageMode CoverageV, typename Binding,
          typename Family = void, ResourceKey Key = "value">
struct CoveragePolicyFor;
/**
 * @brief Full coverage.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 * @tparam Key Value resource key.
 */
template <typename B, typename Family, ResourceKey Key>
struct CoveragePolicyFor<ProjectionCoverageMode::NONE, B, Family, Key> {
  using Type = Pullback::ProjectionCoverage::None;
};
/**
 * @brief Coverage equal to the projection value weight.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 * @tparam Key Value resource key.
 */
template <typename B, typename Family, ResourceKey Key>
struct CoveragePolicyFor<ProjectionCoverageMode::WEIGHT, B, Family, Key> {
  using Type = Pullback::ProjectionCoverage::Weight;
};
/**
 * @brief Coverage equal to the squared value weight.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 * @tparam Key Value resource key.
 */
template <typename B, typename Family, ResourceKey Key>
struct CoveragePolicyFor<ProjectionCoverageMode::WEIGHT_SQUARED, B, Family,
                         Key> {
  using Type = Pullback::ProjectionCoverage::WeightSquared;
};
/**
 * @brief Coverage faded across the projection edge.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 * @tparam Key Value resource key.
 */
template <typename B, typename Family, ResourceKey Key>
struct CoveragePolicyFor<ProjectionCoverageMode::EDGE_FADE, B, Family, Key> {
  using Type =
      Pullback::ProjectionCoverage::EdgeFade<ValueProvider<B, Family, Key>>;
};

/** @brief Squared projection edge-fade coverage.
 * @tparam B Frame binding.
 * @tparam Family Edge-width family.
 * @tparam Key Parameter resource key. */
template <typename B, typename Family, ResourceKey Key>
struct CoveragePolicyFor<ProjectionCoverageMode::EDGE_FADE_SQUARED, B, Family,
                         Key> {
  using Type =
      ProjectionCoverage::EdgeFadeSquared<ValueProvider<B, Family, Key>>;
};
