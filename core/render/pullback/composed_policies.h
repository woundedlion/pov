/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/**
 * @brief Metadata accompanying an effect's explicit ranked stage pipeline.
 * @details A derived Spec supplies `template <typename B> using Pipeline`.
 * PROJECTION controls projection sliders; TRANSFER, COVERAGE and FIELD_COVERAGE
 * describe material stages. LensPolicy and SURFACE_PLACEMENT describe lens and
 * displacement ordering. HARMONY, HUE and BRIGHTNESS select color behavior;
 * ANIMATED_PROJECTION controls projection clocks. The *PolicyFor helpers are
 * optional conveniences for authoring the Pipeline alias.
 */
struct Spec {
  static constexpr ProjectionKind PROJECTION = ProjectionKind::STEREOGRAPHIC;
  static constexpr TransferKind TRANSFER = TransferKind::NONE;
  static constexpr ProjectionCoverageMode COVERAGE =
      ProjectionCoverageMode::WEIGHT;
  static constexpr FieldCoverageKind FIELD_COVERAGE = FieldCoverageKind::NONE;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  static constexpr HueMode HUE = HueMode::NONE;
  static constexpr Color::BrightnessEnvelope BRIGHTNESS =
      Color::BrightnessEnvelope::NONE;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr SurfacePlacement SURFACE_PLACEMENT =
      SurfacePlacement::BEFORE_LENS;
  using LensPolicy = void;
};

template <typename Family, typename Binding> struct SourcePolicyFor;
template <typename B> struct SourcePolicyFor<GridSourceParams, B> {
  using Type = Pullback::Source::Grid<SourceProvider<B, GridSourceParams>>;
};
template <typename B> struct SourcePolicyFor<TwinWaveSourceParams, B> {
  using Type =
      Pullback::Source::TwinWave<SourceProvider<B, TwinWaveSourceParams>>;
};
template <typename B> struct SourcePolicyFor<SpiralSourceParams, B> {
  using Type = Pullback::Source::Spiral<SourceProvider<B, SpiralSourceParams>>;
};
template <typename B> struct SourcePolicyFor<LatticeSourceParams, B> {
  using Type = Pullback::Source::PrimitiveLattice<
      SourceProvider<B, LatticeSourceParams>>;
};
template <typename B> struct SourcePolicyFor<ProjectedNoiseSourceParams, B> {
  using Type = Pullback::Source::ProjectedNoise<
      SourceProvider<B, ProjectedNoiseSourceParams>, math::NoiseBasis::SIMPLEX>;
};
template <typename B> struct SourcePolicyFor<SphericalNoiseSourceParams, B> {
  using Type = Pullback::Source::SphericalNoise<
      SourceProvider<B, SphericalNoiseSourceParams>, math::NoiseBasis::SIMPLEX>;
};

template <typename Family, typename Binding, ResourceKey Key, bool TrackPath>
struct WarpPolicyFor;
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<MirrorParams, B, K, T> {
  using Type = Pullback::Warp::MirrorTile<WarpProvider<B, K, MirrorParams, T>>;
};
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<WaveShearParams, B, K, T> {
  using Type =
      Pullback::Warp::WaveShear<WarpProvider<B, K, WaveShearParams, T>>;
};
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<VectorNoiseParams, B, K, T> {
  using Type =
      Pullback::Warp::VectorNoise<WarpProvider<B, K, VectorNoiseParams, T>,
                                  math::NoiseBasis::SIMPLEX,
                                  Pullback::Warp::FlatEnvelope>;
};
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<AffineParams, B, K, T> {
  using Type = Pullback::Warp::AffineFrame<WarpProvider<B, K, AffineParams, T>>;
};
template <typename B, ResourceKey K, bool T>
struct WarpPolicyFor<PolarParams, B, K, T> {
  using Type = Pullback::Warp::PolarChart<WarpProvider<B, K, PolarParams, T>,
                                          Pullback::Warp::LinearPolar, 1>;
};

template <typename Family, typename Binding, bool TrackPath>
struct SurfacePolicyFor;
template <typename B, bool T>
struct SurfacePolicyFor<SurfaceNoiseParams, B, T> {
  using Type =
      Pullback::Surface::CurlNoise<SurfaceProvider<B, SurfaceNoiseParams, T>,
                                   math::NoiseBasis::SIMPLEX,
                                   Pullback::Surface::Euler>;
};
template <typename B, bool T>
struct SurfacePolicyFor<DirectSurfaceParams, B, T> {
  using Type =
      Pullback::Surface::DirectNoise<SurfaceProvider<B, DirectSurfaceParams, T>,
                                     math::NoiseBasis::SIMPLEX>;
};
template <typename B, bool T>
struct SurfacePolicyFor<PeriodicRippleParams, B, T> {
  using Type = Pullback::Surface::PeriodicRipple<
      SurfaceProvider<B, PeriodicRippleParams, T>>;
};

template <ProjectionKind ProjectionV, typename Binding>
struct ProjectionPolicyFor;
template <typename B>
struct ProjectionPolicyFor<ProjectionKind::STEREOGRAPHIC, B> {
  using Type = Pullback::Projection::Stereographic<ProjectionProvider<B>>;
};
template <typename B>
struct ProjectionPolicyFor<ProjectionKind::GNOMONIC_FOLDED, B> {
  using Type = Pullback::Projection::Gnomonic<
      ProjectionProvider<B>, Pullback::Projection::GnomonicHemisphere::FOLDED>;
};
template <typename B>
struct ProjectionPolicyFor<ProjectionKind::EQUIRECTANGULAR, B> {
  using Type = Pullback::Projection::Equirectangular<ProjectionProvider<B>>;
};
template <typename B>
struct ProjectionPolicyFor<ProjectionKind::FOLDED_SINUSOIDAL, B> {
  using Type = Pullback::Projection::FoldedSinusoidal<ProjectionProvider<B>>;
};

template <TransferKind TransferV, typename Binding, typename Family = void>
struct TransferPolicyFor;
template <typename B, typename Family>
struct TransferPolicyFor<TransferKind::ISO_CONTOUR, B, Family> {
  using Type = Pullback::Transfer::IsoContour<ValueProvider<B, Family>>;
};

template <TransferKind TransferV, typename Binding, typename Family = void>
struct TransferStageFor;
template <typename B, typename Family>
struct TransferStageFor<TransferKind::NONE, B, Family> {
  using Type = void;
};
template <typename B, typename Family>
struct TransferStageFor<TransferKind::ISO_CONTOUR, B, Family> {
  using Type = Pullback::Stage::Transfer<
      typename TransferPolicyFor<TransferKind::ISO_CONTOUR, B, Family>::Type>;
};

template <ProjectionCoverageMode CoverageV, typename Binding,
          typename Family = void>
struct CoveragePolicyFor;
template <typename B, typename Family>
struct CoveragePolicyFor<ProjectionCoverageMode::NONE, B, Family> {
  using Type = Pullback::ProjectionCoverage::None;
};
template <typename B, typename Family>
struct CoveragePolicyFor<ProjectionCoverageMode::WEIGHT, B, Family> {
  using Type = Pullback::ProjectionCoverage::Weight;
};
template <typename B, typename Family>
struct CoveragePolicyFor<ProjectionCoverageMode::WEIGHT_SQUARED, B, Family> {
  using Type = Pullback::ProjectionCoverage::WeightSquared;
};
template <typename B, typename Family>
struct CoveragePolicyFor<ProjectionCoverageMode::EDGE_FADE, B, Family> {
  using Type = Pullback::ProjectionCoverage::EdgeFade<ValueProvider<B, Family>>;
};

template <FieldCoverageKind CoverageV, typename Binding, typename Family = void>
struct FieldCoverageStageFor;
template <typename B, typename Family>
struct FieldCoverageStageFor<FieldCoverageKind::NONE, B, Family> {
  using Type = void;
};
template <typename B, typename Family>
struct FieldCoverageStageFor<FieldCoverageKind::VALUE_CUTOUT, B, Family> {
  using Type = Pullback::Stage::ApplyCoverage<
      Pullback::ValueCoverage::ValueCutout<ValueProvider<B, Family>>>;
};
