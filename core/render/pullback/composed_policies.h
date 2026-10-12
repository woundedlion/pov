/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/** @file composed_policies.h
 * @brief Composed-effect `Spec` metadata and the policy selectors that map
 * it to pipeline stages.
 */

/** @brief Authored controls accompanying an explicit typed stage pipeline.
 * @details Derived specs supply `template <typename B> using Pipeline`.
 * Stage topology is exposed by `MetadataFor`. */
struct Spec {
  /// Hue harmony of generated palettes.
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  /// Whether the projection frame random-walks and spins.
  static constexpr bool ANIMATED_PROJECTION = true;
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
