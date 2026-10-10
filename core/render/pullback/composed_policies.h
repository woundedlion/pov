/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/** @file composed_policies.h
 * @brief Composed-effect `Spec` metadata and the policy selectors that map
 * it to pipeline stages.
 */

/**
 * @brief Metadata accompanying an effect's explicit ranked stage pipeline.
 * @details A derived Spec supplies `template <typename B> using Pipeline`;
 * the *PolicyFor helpers are optional conveniences for authoring it.
 */
struct Spec {
  /**
   * @brief Sphere-to-plane projection of the pipeline's `Stage::Project`.
   * @details Must name the policy `ProjectionPolicyFor<PROJECTION, B>::Type`
   * the pipeline projects through; a pipeline without a Project stage fails
   * to compile. Also gates the central-meridian (`uses_central_meridian`) and
   * singularity-fade (`uses_singularity_fade`) projection parameters.
   * Default: stereographic.
   */
  static constexpr ProjectionKind PROJECTION = ProjectionKind::STEREOGRAPHIC;
  /**
   * @brief Transfer curve the pipeline applies to the sampled value.
   * @details `NONE` requires a pipeline with no `Stage::Transfer`; any other
   * kind requires one. `ISO_CONTOUR` also requires `iso_level` and `iso_width`
   * in the "value" parameter family. Default: no transfer stage.
   */
  static constexpr TransferKind TRANSFER = TransferKind::NONE;
  /**
   * @brief Projection coverage the `Stage::Sample` crossing applies.
   * @details Must match the Sample stage's coverage policy (`NONE` when the
   * pipeline has no Sample stage). `EDGE_FADE` also requires an `edge_width`
   * "value" field and a projection that reports edge distance. Default:
   * coverage equals the projection's value weight.
   */
  static constexpr ProjectionCoverageMode COVERAGE =
      ProjectionCoverageMode::WEIGHT;
  /**
   * @brief Value-dependent coverage applied after sampling.
   * @details `NONE` requires a pipeline with no `Stage::ApplyCoverage`; any
   * other kind requires one. `VALUE_CUTOUT` also requires `cutout_threshold`
   * and `cutout_softness` in the "value" family. Default: no cutout stage.
   */
  static constexpr FieldCoverageKind FIELD_COVERAGE = FieldCoverageKind::NONE;
  /** @brief Hue harmony of every palette in the generated palette cycle.
      Default: triadic. */
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  /**
   * @brief What drives the colour stage's per-pixel hue rotation.
   * @details Gates the hue-rotation colour parameters and the hue LUT storage
   * (`NOISE` adds the hue-noise field and its LUT). `PATH_LENGTH` must be set
   * exactly when a pipeline warp or surface provider tracks path length.
   * Default: no hue rotation.
   */
  static constexpr HueMode HUE = HueMode::NONE;
  /**
   * @brief Brightness envelope the colour stage applies over the palette
   *        coordinate.
   * @details Gates the brightness-range colour parameters; any value other than
   * `NONE` exposes them. Default: flat, full brightness.
   */
  static constexpr Color::BrightnessEnvelope BRIGHTNESS =
      Color::BrightnessEnvelope::NONE;
  /**
   * @brief Whether the projection frame random-walks and spins.
   * @details True adds the projection walk, its noise and the spin-rate and
   * wander projection parameters; false pins the projection to the identity
   * frame. Default: animated.
   */
  static constexpr bool ANIMATED_PROJECTION = true;
  /**
   * @brief Pipeline order of the surface (`Stage::Displace`) stage relative to
   *        the lens stage.
   * @details Checked only when the pipeline has both stages. Default: surface
   * displacement before the lens.
   */
  static constexpr SurfacePlacement SURFACE_PLACEMENT =
      SurfacePlacement::BEFORE_LENS;
  /**
   * @brief Policy type of the pipeline's non-Mobius `Stage::Lens`.
   * @details Must equal that stage's `LensPolicy`; `void` when the pipeline has
   * no lens stage or only a `Lens::Mobius` one. Default: `void`.
   */
  using LensPolicy = void;
};

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
/** @brief Spiral source. @tparam B Frame binding. */
template <typename B> struct SourcePolicyFor<SpiralSourceParams, B> {
  using Type = Pullback::Source::Spiral<SourceProvider<B, SpiralSourceParams>>;
};
/** @brief Primitive-lattice source. @tparam B Frame binding. */
template <typename B> struct SourcePolicyFor<LatticeSourceParams, B> {
  using Type = Pullback::Source::PrimitiveLattice<
      SourceProvider<B, LatticeSourceParams>>;
};
/** @brief Planar simplex-noise source. @tparam B Frame binding. */
template <typename B> struct SourcePolicyFor<ProjectedNoiseSourceParams, B> {
  using Type = Pullback::Source::ProjectedNoise<
      SourceProvider<B, ProjectedNoiseSourceParams>, math::NoiseBasis::SIMPLEX>;
};
/** @brief Spherical simplex-noise source. @tparam B Frame binding. */
template <typename B> struct SourcePolicyFor<SphericalNoiseSourceParams, B> {
  using Type = Pullback::Source::SphericalNoise<
      SourceProvider<B, SphericalNoiseSourceParams>, math::NoiseBasis::SIMPLEX>;
};

/**
 * @brief Maps a warp parameter family to its `Pullback::Warp` policy, exposed
 *        as `Type`.
 * @tparam Family Warp parameter family.
 * @tparam Binding Frame binding of the effect.
 * @tparam Key Parameter resource key of this warp instance.
 * @tparam TrackPath Whether the warp accumulates path length for hue.
 */
template <typename Family, typename Binding, ResourceKey Key, bool TrackPath>
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
 */
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<AffineParams, B, K, T> {
  using Type = Pullback::Warp::AffineFrame<WarpProvider<B, K, AffineParams, T>>;
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

/**
 * @brief Maps a `TransferKind` to its `Pullback::Transfer` policy, exposed as
 *        `Type`; undefined for `TransferKind::NONE`.
 * @tparam TransferV Transfer kind.
 * @tparam Binding Frame binding of the effect.
 * @tparam Family Value parameter family; void when none is read.
 */
template <TransferKind TransferV, typename Binding, typename Family = void>
struct TransferPolicyFor;
/**
 * @brief Iso-contour transfer.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 */
template <typename B, typename Family>
struct TransferPolicyFor<TransferKind::ISO_CONTOUR, B, Family> {
  using Type = Pullback::Transfer::IsoContour<ValueProvider<B, Family>>;
};

/**
 * @brief Maps a `TransferKind` to its `Stage::Transfer`, exposed as `Type`
 *        (void for `TransferKind::NONE`).
 * @tparam TransferV Transfer kind.
 * @tparam Binding Frame binding of the effect.
 * @tparam Family Value parameter family; void when none is read.
 */
template <TransferKind TransferV, typename Binding, typename Family = void>
struct TransferStageFor;
/**
 * @brief No transfer stage.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 */
template <typename B, typename Family>
struct TransferStageFor<TransferKind::NONE, B, Family> {
  using Type = void;
};
/**
 * @brief Iso-contour transfer stage.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 */
template <typename B, typename Family>
struct TransferStageFor<TransferKind::ISO_CONTOUR, B, Family> {
  using Type = Pullback::Stage::Transfer<
      typename TransferPolicyFor<TransferKind::ISO_CONTOUR, B, Family>::Type>;
};

/**
 * @brief Maps a `ProjectionCoverageMode` to its `Pullback::ProjectionCoverage`
 *        policy, exposed as `Type`.
 * @tparam CoverageV Coverage mode.
 * @tparam Binding Frame binding of the effect.
 * @tparam Family Value parameter family; void when none is read.
 */
template <ProjectionCoverageMode CoverageV, typename Binding,
          typename Family = void>
struct CoveragePolicyFor;
/**
 * @brief Full coverage.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 */
template <typename B, typename Family>
struct CoveragePolicyFor<ProjectionCoverageMode::NONE, B, Family> {
  using Type = Pullback::ProjectionCoverage::None;
};
/**
 * @brief Coverage equal to the projection value weight.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 */
template <typename B, typename Family>
struct CoveragePolicyFor<ProjectionCoverageMode::WEIGHT, B, Family> {
  using Type = Pullback::ProjectionCoverage::Weight;
};
/**
 * @brief Coverage equal to the squared value weight.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 */
template <typename B, typename Family>
struct CoveragePolicyFor<ProjectionCoverageMode::WEIGHT_SQUARED, B, Family> {
  using Type = Pullback::ProjectionCoverage::WeightSquared;
};
/**
 * @brief Coverage faded across the projection edge.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 */
template <typename B, typename Family>
struct CoveragePolicyFor<ProjectionCoverageMode::EDGE_FADE, B, Family> {
  using Type = Pullback::ProjectionCoverage::EdgeFade<ValueProvider<B, Family>>;
};

/**
 * @brief Maps a `FieldCoverageKind` to its `Stage::ApplyCoverage`, exposed as
 *        `Type` (void for `FieldCoverageKind::NONE`).
 * @tparam CoverageV Field coverage kind.
 * @tparam Binding Frame binding of the effect.
 * @tparam Family Value parameter family; void when none is read.
 */
template <FieldCoverageKind CoverageV, typename Binding, typename Family = void>
struct FieldCoverageStageFor;
/**
 * @brief No field coverage stage.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 */
template <typename B, typename Family>
struct FieldCoverageStageFor<FieldCoverageKind::NONE, B, Family> {
  using Type = void;
};
/**
 * @brief Value-cutout coverage stage.
 * @tparam B Frame binding.
 * @tparam Family Value parameter family.
 */
template <typename B, typename Family>
struct FieldCoverageStageFor<FieldCoverageKind::VALUE_CUTOUT, B, Family> {
  using Type = Pullback::Stage::ApplyCoverage<
      Pullback::ValueCoverage::ValueCutout<ValueProvider<B, Family>>>;
};
