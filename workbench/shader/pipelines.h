/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_SHADER_WORKBENCH

/**
 * @file pipelines.h
 * @brief The compiled inverse-pipeline catalog: the stage adapters each slot tuple selects, the pipelines they compose into, and the manifest that matches a configuration to one.
 */

#include "workbench/shader/bindings.h"
#include "workbench/shader/presets.h"

namespace Workbench {

using CodeEmission = Pullback::CodeEmission;
using ApproximationOracleId = Pullback::ApproximationOracleId;
using ApproximationDomain = Pullback::ApproximationDomain;
using ApproximationAggregation = Pullback::ApproximationAggregation;
using ApproximationMetric = Pullback::ApproximationMetric;

struct TopologyKey {
  Function function{};
  Projection projection{};
  ProjectionFramePolicy projection_frame{};
  SurfaceLens surface_lens{};
  SignalWeight signal_weight{};
  ValueTransfer value_transfer{};
  CoveragePolicy coverage{};
  PeirceLayout peirce_layout{};
  AiroceanLayout airocean_layout{};
  BonneHemisphere bonne_hemisphere{};
  GnomonicHemispherePolicy gnomonic_hemisphere{};
  SurfaceNoise surface_noise{};
  SurfaceNoisePlacement surface_noise_placement{};
  NoiseBasis surface_noise_basis{};
  SurfaceCurlIntegrator surface_curl_integrator{};
  NoiseBasis source_noise_basis{};
  WarpStageKind outer_warp{};
  NoiseBasis outer_warp_basis{};
  WarpEnvelope outer_warp_envelope{};
  PolarMode outer_polar_mode{};
  CurlIntegrator outer_curl_integrator{};
  uint8_t outer_polar_harmonic{};
  WarpStageKind inner_warp{};
  NoiseBasis inner_warp_basis{};
  WarpEnvelope inner_warp_envelope{};
  PolarMode inner_polar_mode{};
  CurlIntegrator inner_curl_integrator{};
  uint8_t inner_polar_harmonic{};

  constexpr bool operator==(const TopologyKey &) const = default;
};

template <typename... Stages>
using InversePipeline = Pullback::Pipeline<ShaderWorkbenchBinding, Stages...>;

#if defined(__IMXRT1062__)
HS_FLASH_MEMBER
#else
__attribute__((always_inline))
#endif
inline Vector pullback_outer_camera_lookup(const Vector &input,
                                           const FrameState &frame) {
  return rotate(input, frame.transforms.outer_conj);
}

struct OuterCameraStage
    : Pullback::Stage::Contract<OuterCameraStage, Pullback::SphereSample,
                                Pullback::SphereSample> {
  using Policies = std::tuple<>;

  static constexpr bool implements(const TopologyKey &) { return true; }

  template <typename Binding>
  __attribute__((always_inline)) static Pullback::SphereSample
  run(const Pullback::SphereSample &input,
      const typename Binding::FrameState &frame, const Pullback::NoPrepared &) {
    return {pullback_outer_camera_lookup(input.dir, frame), input.path_length};
  }
};

template <SurfaceLens LensV>
using LensPolicy = std::conditional_t<
    LensV == SurfaceLens::NONE, void,
    std::conditional_t<
        LensV == SurfaceLens::GLITCH, Pullback::Lens::Glitch,
        std::conditional_t<
            LensV == SurfaceLens::KALEIDOSCOPE, Pullback::Lens::Kaleidoscope,
            std::conditional_t<
                LensV == SurfaceLens::MOBIUS,
                Pullback::Lens::Mobius<LensStateProvider>,
                std::conditional_t<
                    LensV == SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
                    Pullback::Lens::DodecahedralKaleidoscope,
                    std::conditional_t<
                        LensV == SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM,
                        Pullback::Lens::HexagonalPrismKaleidoscope,
                        std::conditional_t<
                            LensV == SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM,
                            Pullback::Lens::PentagonalPrismKaleidoscope,
                            Pullback::Lens::TriangularPrismKaleidoscope>>>>>>>;

template <Projection ProjectionV>
using ProjectionPolicy = std::conditional_t<
    ProjectionV == Projection::STEREOGRAPHIC,
    Pullback::Projection::Stereographic<ProjectionStateProvider>,
    std::conditional_t<
        ProjectionV == Projection::GNOMONIC,
        Pullback::Projection::Gnomonic<
            ProjectionStateProvider,
            Pullback::Projection::GnomonicHemisphere::FOLDED>,
        std::conditional_t<
            ProjectionV == Projection::BONNE,
            Pullback::Projection::Bonne<ProjectionStateProvider, true>,
            std::conditional_t<
                ProjectionV == Projection::EQUIRECTANGULAR,
                Pullback::Projection::Equirectangular<ProjectionStateProvider>,
                Pullback::Projection::PeirceSquare<ProjectionStateProvider>>>>>;

template <SurfaceLens LensV>
struct SelectedLensStage : Pullback::Stage::Lens<LensPolicy<LensV>> {
  static_assert(LensV == SurfaceLens::GLITCH ||
                LensV == SurfaceLens::KALEIDOSCOPE ||
                LensV == SurfaceLens::MOBIUS ||
                LensV == SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL ||
                LensV == SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM ||
                LensV == SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM ||
                LensV == SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM);
  static constexpr bool implements(const TopologyKey &key) {
    return key.surface_lens == LensV;
  }
};

template <Projection ProjectionV, SurfaceLens LensV>
struct SelectedProjectStage
    : Pullback::Stage::Project<ProjectionPolicy<ProjectionV>> {
  static_assert(ProjectionV == Projection::STEREOGRAPHIC ||
                ProjectionV == Projection::GNOMONIC ||
                ProjectionV == Projection::BONNE ||
                ProjectionV == Projection::EQUIRECTANGULAR ||
                ProjectionV == Projection::PEIRCE_QUINCUNCIAL);
  static constexpr bool implements(const TopologyKey &key) {
    return key.projection == ProjectionV && key.surface_lens == LensV &&
           key.surface_noise == SurfaceNoise::NONE &&
           (ProjectionV != Projection::GNOMONIC ||
            key.gnomonic_hemisphere == GnomonicHemispherePolicy::FOLDED) &&
           (ProjectionV != Projection::PEIRCE_QUINCUNCIAL ||
            key.peirce_layout == PeirceLayout::SQUARE) &&
           (ProjectionV != Projection::BONNE ||
            key.bonne_hemisphere == BonneHemisphere::NORTH);
  }
};

struct SinusoidalCurlDisplaceStage
    : Pullback::Stage::Displace<Pullback::Surface::CurlNoise<
          SurfaceStateProvider, NoiseBasis::SIMPLEX,
          Pullback::Surface::Euler>> {
  static constexpr bool implements(const TopologyKey &key) {
    return key.surface_noise == SurfaceNoise::CURL &&
           key.surface_noise_placement == SurfaceNoisePlacement::BEFORE_LENS &&
           key.surface_noise_basis == NoiseBasis::SIMPLEX &&
           key.surface_curl_integrator == SurfaceCurlIntegrator::EULER;
  }
};

struct SinusoidalCurlProjectStage
    : Pullback::Stage::Project<
          Pullback::Projection::FoldedSinusoidal<ProjectionStateProvider>> {
  static constexpr bool implements(const TopologyKey &key) {
    return key.projection == Projection::SINUSOIDAL &&
           key.surface_lens == SurfaceLens::NONE;
  }
};

using SinusoidalCurlSphereRun =
    Pullback::Stage::Placed<Pullback::CodeEmission::OUT_OF_LINE_FLASH,
                            SinusoidalCurlDisplaceStage,
                            SinusoidalCurlProjectStage>;

template <WarpStageKind KindV, bool Outer>
using WarpPolicy = std::conditional_t<
    KindV == WarpStageKind::NONE, void,
    std::conditional_t<
        KindV == WarpStageKind::AFFINE_FRAME,
        Pullback::Warp::AffineFrame<WarpStateProvider<Outer>>,
        std::conditional_t<
            KindV == WarpStageKind::WAVE_SHEAR,
            Pullback::Warp::WaveShear<WarpStateProvider<Outer>>,
            std::conditional_t<
                KindV == WarpStageKind::VECTOR_NOISE,
                Pullback::Warp::VectorNoise<WarpStateProvider<Outer>,
                                            NoiseBasis::SIMPLEX,
                                            Pullback::Warp::FlatEnvelope>,
                std::conditional_t<
                    KindV == WarpStageKind::MIRROR_TILE,
                    Pullback::Warp::MirrorTile<WarpStateProvider<Outer>>,
                    Pullback::Warp::PolarChart<WarpStateProvider<Outer>,
                                               Pullback::Warp::LinearPolar,
                                               1>>>>>>;

template <WarpStageKind KindV, bool Outer>
struct SelectedWarpStage : Pullback::Stage::Warp<WarpPolicy<KindV, Outer>> {
  static_assert(KindV == WarpStageKind::AFFINE_FRAME ||
                KindV == WarpStageKind::WAVE_SHEAR ||
                KindV == WarpStageKind::VECTOR_NOISE ||
                KindV == WarpStageKind::MIRROR_TILE ||
                KindV == WarpStageKind::POLAR_CHART);
  static_assert(Outer || KindV == WarpStageKind::WAVE_SHEAR ||
                KindV == WarpStageKind::MIRROR_TILE);
  static constexpr bool implements(const TopologyKey &key) {
    if constexpr (Outer)
      return key.outer_warp == KindV &&
             (KindV != WarpStageKind::WAVE_SHEAR ||
              key.outer_warp_envelope == WarpEnvelope::FLAT);
    else
      return key.inner_warp == KindV &&
             (KindV != WarpStageKind::WAVE_SHEAR ||
              key.inner_warp_envelope == WarpEnvelope::FLAT);
  }
};

template <Function FunctionV>
using SourcePolicy = std::conditional_t<
    FunctionV == Function::GRID, Pullback::Source::Grid<SourceStateProvider>,
    std::conditional_t<FunctionV == Function::PRIMITIVE_LATTICE,
                       Pullback::Source::PrimitiveLattice<SourceStateProvider>,
                       Pullback::Source::TwinWave<SourceStateProvider>>>;

template <CoveragePolicy CoverageV> struct ProjectionCoverageMapping {
  static_assert(CoverageV == CoveragePolicy::OPAQUE ||
                CoverageV == CoveragePolicy::PROJECTION_WEIGHT ||
                CoverageV == CoveragePolicy::PROJECTION_WEIGHT_SQUARED ||
                CoverageV == CoveragePolicy::EDGE_FADE);
  static constexpr Pullback::ProjectionCoverageMode MODE =
      CoverageV == CoveragePolicy::OPAQUE
          ? Pullback::ProjectionCoverageMode::NONE
      : CoverageV == CoveragePolicy::PROJECTION_WEIGHT
          ? Pullback::ProjectionCoverageMode::WEIGHT
      : CoverageV == CoveragePolicy::PROJECTION_WEIGHT_SQUARED
          ? Pullback::ProjectionCoverageMode::WEIGHT_SQUARED
          : Pullback::ProjectionCoverageMode::EDGE_FADE;
  using Type = std::conditional_t<
      MODE == Pullback::ProjectionCoverageMode::NONE,
      Pullback::ProjectionCoverage::None,
      std::conditional_t<
          MODE == Pullback::ProjectionCoverageMode::WEIGHT,
          Pullback::ProjectionCoverage::Weight,
          std::conditional_t<
              MODE == Pullback::ProjectionCoverageMode::WEIGHT_SQUARED,
              Pullback::ProjectionCoverage::WeightSquared,
              Pullback::ProjectionCoverage::EdgeFade<ValueStateProvider>>>>;
};

template <CoveragePolicy CoverageV>
using CoveragePolicyFor = typename ProjectionCoverageMapping<CoverageV>::Type;

template <Function FunctionV, ValueTransfer TransferV, CoveragePolicy CoverageV>
struct SelectedSampleStage
    : Pullback::Stage::Sample<SourcePolicy<FunctionV>,
                              Pullback::Weight::Projection,
                              CoveragePolicyFor<CoverageV>> {
  static_assert(FunctionV == Function::GRID ||
                FunctionV == Function::PRIMITIVE_LATTICE ||
                FunctionV == Function::TWIN_WAVE);
  static_assert(CoverageV == CoveragePolicy::OPAQUE ||
                CoverageV == CoveragePolicy::EDGE_FADE ||
                CoverageV == CoveragePolicy::PROJECTION_WEIGHT_SQUARED ||
                CoverageV == CoveragePolicy::PROJECTION_WEIGHT);
  static constexpr CoveragePolicy COVERAGE = CoverageV;
  static constexpr bool implements(const TopologyKey &key) {
    return key.function == FunctionV &&
           key.signal_weight == SignalWeight::PROJECTION &&
           key.value_transfer == TransferV && key.coverage == CoverageV;
  }
};

struct IsoContourTransferStage
    : Pullback::Stage::Transfer<
          Pullback::Transfer::IsoContour<ValueStateProvider>> {
  static constexpr bool implements(const TopologyKey &key) {
    return key.value_transfer == ValueTransfer::ISO_CONTOUR;
  }
};

struct ColorStage : Pullback::Stage::Colorize<
                        Pullback::Color::GeneratedPalette<ColorStateProvider>> {
  static constexpr bool implements(const TopologyKey &) { return true; }
};

using GlitchNoiseGridWaveShearPipelineBase = InversePipeline<
    OuterCameraStage, SelectedLensStage<SurfaceLens::GLITCH>,
    SelectedProjectStage<Projection::STEREOGRAPHIC, SurfaceLens::GLITCH>,
    SelectedWarpStage<WarpStageKind::WAVE_SHEAR, true>,
    SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                        CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
    ColorStage>;
struct GlitchNoiseGridWaveShearPipeline : GlitchNoiseGridWaveShearPipelineBase {
  HS_HOT_FLASH_MEMBER static Color4 shade_prepared(const Vector &view,
                                                   const FrameState &frame,
                                                   const void *storage) {
    return GlitchNoiseGridWaveShearPipelineBase::evaluate(
        view, frame,
        *static_cast<const typename GlitchNoiseGridWaveShearPipelineBase::
                         PreparedTuple *>(storage));
  }
};
using KaleidoscopeTwinWaveInnerMirrorPipeline = InversePipeline<
    OuterCameraStage, SelectedLensStage<SurfaceLens::KALEIDOSCOPE>,
    SelectedProjectStage<Projection::STEREOGRAPHIC, SurfaceLens::KALEIDOSCOPE>,
    SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
    SelectedSampleStage<Function::TWIN_WAVE, ValueTransfer::NONE,
                        CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
    ColorStage>;
using StereographicMobiusTwinWaveInnerMirrorPipeline = InversePipeline<
    OuterCameraStage, SelectedLensStage<SurfaceLens::MOBIUS>,
    SelectedProjectStage<Projection::STEREOGRAPHIC, SurfaceLens::MOBIUS>,
    SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
    SelectedSampleStage<Function::TWIN_WAVE, ValueTransfer::NONE,
                        CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
    ColorStage>;
using StereographicHexagonalPrismTwinWaveInnerMirrorPipeline = InversePipeline<
    OuterCameraStage,
    SelectedLensStage<SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM>,
    SelectedProjectStage<Projection::STEREOGRAPHIC,
                         SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM>,
    SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
    SelectedSampleStage<Function::TWIN_WAVE, ValueTransfer::NONE,
                        CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
    ColorStage>;
using GnomonicKaleidoscopeGridMirrorPipeline = InversePipeline<
    OuterCameraStage, SelectedLensStage<SurfaceLens::KALEIDOSCOPE>,
    SelectedProjectStage<Projection::GNOMONIC, SurfaceLens::KALEIDOSCOPE>,
    SelectedWarpStage<WarpStageKind::MIRROR_TILE, true>,
    SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                        CoveragePolicy::EDGE_FADE>,
    ColorStage>;
using GnomonicAlienCoreMirrorPipeline = InversePipeline<
    OuterCameraStage, SelectedLensStage<SurfaceLens::GLITCH>,
    SelectedProjectStage<Projection::GNOMONIC, SurfaceLens::GLITCH>,
    SelectedWarpStage<WarpStageKind::MIRROR_TILE, true>,
    SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                        CoveragePolicy::EDGE_FADE>,
    ColorStage>;
using PeirceDodecahedralGridPipelineBase = InversePipeline<
    OuterCameraStage, SelectedLensStage<SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
    SelectedProjectStage<Projection::PEIRCE_QUINCUNCIAL,
                         SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
    SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                        CoveragePolicy::EDGE_FADE>,
    ColorStage>;
struct PeirceDodecahedralGridPipeline : PeirceDodecahedralGridPipelineBase {
  HS_HOT_FLASH_MEMBER static Color4 shade_prepared(const Vector &view,
                                                   const FrameState &frame,
                                                   const void *storage) {
    return PeirceDodecahedralGridPipelineBase::evaluate(
        view, frame,
        *static_cast<
            const typename PeirceDodecahedralGridPipelineBase::PreparedTuple *>(
            storage));
  }
};
using GnomonicDodecahedralGridWaveMirrorPipelineBase = InversePipeline<
    OuterCameraStage, SelectedLensStage<SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
    SelectedProjectStage<Projection::GNOMONIC,
                         SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
    SelectedWarpStage<WarpStageKind::WAVE_SHEAR, true>,
    SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
    SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                        CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
    ColorStage>;
struct GnomonicDodecahedralGridWaveMirrorPipeline
    : GnomonicDodecahedralGridWaveMirrorPipelineBase {
  HS_HOT_FLASH_MEMBER static Color4 shade_prepared(const Vector &view,
                                                   const FrameState &frame,
                                                   const void *storage) {
    return GnomonicDodecahedralGridWaveMirrorPipelineBase::evaluate(
        view, frame,
        *static_cast<
            const typename GnomonicDodecahedralGridWaveMirrorPipelineBase::
                PreparedTuple *>(storage));
  }
};
using GnomonicDodecahedralGridVectorMirrorPipeline = InversePipeline<
    OuterCameraStage, SelectedLensStage<SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
    SelectedProjectStage<Projection::GNOMONIC,
                         SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
    SelectedWarpStage<WarpStageKind::VECTOR_NOISE, true>,
    SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
    SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                        CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
    ColorStage>;
using GnomonicAffineLatticeContourPipeline = InversePipeline<
    OuterCameraStage,
    SelectedProjectStage<Projection::GNOMONIC, SurfaceLens::NONE>,
    SelectedWarpStage<WarpStageKind::AFFINE_FRAME, true>,
    SelectedSampleStage<Function::PRIMITIVE_LATTICE, ValueTransfer::ISO_CONTOUR,
                        CoveragePolicy::PROJECTION_WEIGHT>,
    IsoContourTransferStage, ColorStage>;
using SinusoidalLatticeMeltPipeline = InversePipeline<
    OuterCameraStage, SinusoidalCurlSphereRun,
    SelectedSampleStage<Function::PRIMITIVE_LATTICE, ValueTransfer::NONE,
                        CoveragePolicy::PROJECTION_WEIGHT>,
    ColorStage>;
using StereographicPrismPolarWaveLatticePipeline = InversePipeline<
    OuterCameraStage,
    SelectedLensStage<SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM>,
    SelectedProjectStage<Projection::STEREOGRAPHIC,
                         SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM>,
    SelectedWarpStage<WarpStageKind::POLAR_CHART, true>,
    SelectedWarpStage<WarpStageKind::WAVE_SHEAR, false>,
    SelectedSampleStage<Function::PRIMITIVE_LATTICE, ValueTransfer::NONE,
                        CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
    ColorStage>;
using StereographicDodecahedralGridInnerMirrorPipeline = InversePipeline<
    OuterCameraStage, SelectedLensStage<SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
    SelectedProjectStage<Projection::STEREOGRAPHIC,
                         SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
    SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
    SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                        CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
    ColorStage>;
using EquirectangularDodecahedralGridInnerMirrorPipeline = InversePipeline<
    OuterCameraStage,
    Pullback::Stage::Placed<
        Pullback::CodeEmission::OUT_OF_LINE_FLASH,
        SelectedLensStage<SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
        SelectedProjectStage<Projection::EQUIRECTANGULAR,
                             SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>>,
    SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
    SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                        CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
    ColorStage>;
using StereographicAlienCoreMirrorPipeline = InversePipeline<
    OuterCameraStage, SelectedLensStage<SurfaceLens::GLITCH>,
    SelectedProjectStage<Projection::STEREOGRAPHIC, SurfaceLens::GLITCH>,
    SelectedWarpStage<WarpStageKind::MIRROR_TILE, true>,
    SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                        CoveragePolicy::EDGE_FADE>,
    ColorStage>;

struct ProgramDescriptor {
  InversePipelineId id;
  TopologyKey key;
  ShadeFunction shade;
  void (*prepare)(const FrameState &, void *);
  bool (*continuous_parameters_supported)(const Config &);
  bool (*resources_ready)(const FrameState &);
};

/** Capacity of the per-frame prepared blob a compiled program renders
    from; make_program() pins every program's tuple under it. */
inline constexpr size_t PREPARED_BLOB_BYTES = 256;

inline constexpr void canonicalize_warp_key(WarpStageKind kind,
                                            NoiseBasis &basis,
                                            WarpEnvelope &envelope,
                                            PolarMode &polar_mode,
                                            CurlIntegrator &curl_integrator,
                                            uint8_t &polar_harmonic) {
  if (!warp_uses_noise(kind))
    basis = {};
  if (!warp_uses_envelope(kind))
    envelope = {};
  if (kind != WarpStageKind::CURL_FLOW)
    curl_integrator = {};
  if (kind != WarpStageKind::POLAR_CHART) {
    polar_mode = {};
    polar_harmonic = 0;
  }
}

inline constexpr TopologyKey make_topology_key(const Config &config) {
  const Slots &slots = config.slots;
  TopologyKey key{
      slots.function,
      slots.projection,
      slots.projection_frame,
      slots.surface_lens,
      slots.signal_weight,
      slots.value_transfer,
      slots.coverage,
      slots.peirce_layout,
      slots.airocean_layout,
      slots.bonne_hemisphere,
      slots.gnomonic_hemisphere,
      slots.surface_noise,
      slots.surface_noise_placement,
      config.params.surface_noise.basis,
      config.params.surface_noise.integrator,
      config.params.source.noise_basis,
      slots.warp_program.outer.kind,
      slots.warp_program.outer.basis,
      slots.warp_program.outer.envelope,
      slots.warp_program.outer.polar_mode,
      slots.warp_program.outer.curl_integrator,
      slots.warp_program.outer.polar_harmonic,
      slots.warp_program.inner.kind,
      slots.warp_program.inner.basis,
      slots.warp_program.inner.envelope,
      slots.warp_program.inner.polar_mode,
      slots.warp_program.inner.curl_integrator,
      slots.warp_program.inner.polar_harmonic,
  };
  if (key.projection != Projection::PEIRCE_QUINCUNCIAL)
    key.peirce_layout = {};
  if (key.projection != Projection::AIROCEAN)
    key.airocean_layout = {};
  if (key.projection != Projection::BONNE)
    key.bonne_hemisphere = {};
  if (key.projection != Projection::GNOMONIC)
    key.gnomonic_hemisphere = {};
  if (key.surface_noise == SurfaceNoise::NONE) {
    key.surface_noise_placement = {};
    key.surface_noise_basis = {};
    key.surface_curl_integrator = {};
  }
  if (!source_uses_noise(key.function))
    key.source_noise_basis = {};
  canonicalize_warp_key(key.outer_warp, key.outer_warp_basis,
                        key.outer_warp_envelope, key.outer_polar_mode,
                        key.outer_curl_integrator, key.outer_polar_harmonic);
  canonicalize_warp_key(key.inner_warp, key.inner_warp_basis,
                        key.inner_warp_envelope, key.inner_polar_mode,
                        key.inner_curl_integrator, key.inner_polar_harmonic);
  return key;
}

inline constexpr bool all_continuous_parameters_supported(const Config &) {
  return true;
}

HS_FLASH_MEMBER inline bool pipeline_resources_ready(const FrameState &frame) {
  if (warp_uses_noise(frame.slots.warp_program.outer.kind) &&
      frame.resources.outer_warp_noise == nullptr)
    return false;
  if (warp_uses_noise(frame.slots.warp_program.inner.kind) &&
      frame.resources.inner_warp_noise == nullptr)
    return false;
  if (is_noise_contour(frame.slots.function) &&
      frame.resources.source_noise == nullptr)
    return false;
  if (frame.slots.surface_noise != SurfaceNoise::NONE &&
      frame.resources.surface_noise == nullptr)
    return false;
  if (frame.resources.generated_palette == nullptr)
    return false;
  if (frame.slots.hue_shift == HueShiftMode::NOISE &&
      frame.params.color.hue_shift_amount != 0.0f &&
      frame.resources.color_noise == nullptr)
    return false;
  if (frame.prepared_hue_rotation.active &&
      frame.prepared_hue_rotation.lut == nullptr)
    return false;
  return !frame.prepared_hue_noise.active ||
         frame.prepared_hue_noise.lut != nullptr;
}

/**
 * @brief Builds one manifest row.
 * @tparam Pipeline Compiled inverse pipeline the row selects.
 * @tparam Id Stable identifier for the row.
 * @tparam Key Topology the row is matched on.
 * @param continuous Predicate over the continuous parameters @p Pipeline
 *        serves.
 * @details Rejects at compile time a pipeline whose stages hardcode a
 * topology facet that @p Key does not carry.
 */
template <typename Pipeline, InversePipelineId Id, TopologyKey Key>
inline constexpr ProgramDescriptor
make_program(bool (*continuous)(const Config &)) {
  static_assert(Pipeline::implements(Key),
                "inverse pipeline does not implement its topology key");
  static_assert(sizeof(typename Pipeline::PreparedTuple) <= PREPARED_BLOB_BYTES,
                "prepared blob capacity exceeded");
  return {Id,
          Key,
          &Pipeline::shade_prepared,
          &Pipeline::prepare_into,
          continuous,
          &pipeline_resources_ready};
}

HS_COLD_MEMBER inline const std::array<ProgramDescriptor,
                                       Workbench::INVERSE_PROGRAM_COUNT> &
inverse_programs() {
  static constexpr std::array<ProgramDescriptor,
                              Workbench::INVERSE_PROGRAM_COUNT>
      PROGRAMS{{
          make_program<GlitchNoiseGridWaveShearPipeline,
                       InversePipelineId::GLITCH_NOISE_GRID_WAVE_SHEAR,
                       make_topology_key(
                           Workbench::wave_shear_generated_preset())>(
              &all_continuous_parameters_supported),
          make_program<KaleidoscopeTwinWaveInnerMirrorPipeline,
                       InversePipelineId::KALEIDOSCOPE_TWIN_WAVE_INNER_MIRROR,
                       make_topology_key(
                           Workbench::kaleidoscope_mirror_preset())>(
              &all_continuous_parameters_supported),
          make_program<
              GnomonicKaleidoscopeGridMirrorPipeline,
              InversePipelineId::GNOMONIC_KALEIDOSCOPE_GRID_MIRROR,
              make_topology_key(
                  Workbench::gnomonic_kaleidoscope_grid_mirror_preset())>(
              &all_continuous_parameters_supported),
          make_program<GnomonicAlienCoreMirrorPipeline,
                       InversePipelineId::GNOMONIC_ALIEN_CORE_MIRROR,
                       make_topology_key(Workbench::gnomonic_grid_mirror_preset(
                           SurfaceLens::GLITCH))>(
              &all_continuous_parameters_supported),
          make_program<PeirceDodecahedralGridPipeline,
                       InversePipelineId::PEIRCE_DODECAHEDRAL_GRID,
                       make_topology_key(
                           Workbench::peirce_dodecahedral_generated_preset())>(
              &all_continuous_parameters_supported),
          make_program<
              GnomonicDodecahedralGridWaveMirrorPipeline,
              InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_WAVE_MIRROR,
              make_topology_key(Workbench::gnomonic_wave_shear_grid_preset())>(
              &all_continuous_parameters_supported),
          make_program<
              GnomonicAffineLatticeContourPipeline,
              InversePipelineId::GNOMONIC_AFFINE_LATTICE_CONTOUR,
              make_topology_key(
                  Workbench::gnomonic_affine_lattice_contour_preset())>(
              &all_continuous_parameters_supported),
          make_program<SinusoidalLatticeMeltPipeline,
                       InversePipelineId::SINUSOIDAL_LATTICE_MELT,
                       make_topology_key(
                           Workbench::sinusoidal_lattice_curl_preset(1.0f))>(
              &all_continuous_parameters_supported),
          make_program<
              StereographicPrismPolarWaveLatticePipeline,
              InversePipelineId::STEREOGRAPHIC_PRISM_POLAR_WAVE_LATTICE,
              make_topology_key(
                  Workbench::stereographic_prism_polar_wave_lattice_preset())>(
              &all_continuous_parameters_supported),
          make_program<
              GnomonicDodecahedralGridVectorMirrorPipeline,
              InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_VECTOR_MIRROR,
              make_topology_key(
                  Workbench::
                      gnomonic_dodecahedral_vector_mirror_grid_preset())>(
              &all_continuous_parameters_supported),
          make_program<
              StereographicDodecahedralGridInnerMirrorPipeline,
              InversePipelineId::STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR,
              make_topology_key(
                  Workbench::
                      stereographic_dodecahedral_grid_inner_mirror_preset())>(
              &all_continuous_parameters_supported),
          make_program<
              StereographicHexagonalPrismTwinWaveInnerMirrorPipeline,
              InversePipelineId::
                  STEREOGRAPHIC_HEXAGONAL_PRISM_TWIN_WAVE_INNER_MIRROR,
              make_topology_key(
                  Workbench::
                      stereographic_hexagonal_prism_twin_wave_mirror_preset())>(
              &all_continuous_parameters_supported),
          make_program<
              EquirectangularDodecahedralGridInnerMirrorPipeline,
              InversePipelineId::EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR,
              make_topology_key(
                  Workbench::
                      equirectangular_dodecahedral_double_mapping_grid_inner_mirror_preset())>(
              &all_continuous_parameters_supported),
          make_program<
              StereographicAlienCoreMirrorPipeline,
              InversePipelineId::STEREOGRAPHIC_ALIEN_CORE_MIRROR,
              make_topology_key(
                  Workbench::stereographic_alien_core_mirror_preset())>(
              &all_continuous_parameters_supported),
          make_program<
              StereographicMobiusTwinWaveInnerMirrorPipeline,
              InversePipelineId::STEREOGRAPHIC_MOBIUS_TWIN_WAVE_INNER_MIRROR,
              make_topology_key(
                  Workbench::
                      stereographic_mobius_twin_wave_inner_mirror_preset())>(
              &all_continuous_parameters_supported),
      }};
  return PROGRAMS;
}

HS_COLD_MEMBER inline const ProgramDescriptor *
find_inverse_program(const Config &config) {
  const TopologyKey key = make_topology_key(config);
  for (const ProgramDescriptor &program : inverse_programs())
    if (program.key == key && program.continuous_parameters_supported(config))
      return &program;
  return nullptr;
}

HS_COLD_MEMBER inline const ProgramDescriptor *
get_inverse_program(InversePipelineId id) {
  for (const ProgramDescriptor &program : inverse_programs())
    if (program.id == id)
      return &program;
  return nullptr;
}

HS_COLD_MEMBER inline InversePipelineId
resolve_pipeline_id(const Config &config) {
  const ProgramDescriptor *program = find_inverse_program(config);
  return program == nullptr ? InversePipelineId::NONE : program->id;
}

HS_COLD_MEMBER inline const ProgramDescriptor *
resolve_inverse_program(const FrameState &frame) {
  const Config config{frame.slots, frame.params};
  const ProgramDescriptor *program = find_inverse_program(config);
  if (program == nullptr || !program->resources_ready(frame))
    return nullptr;
  return program;
}

} // namespace Workbench

#endif // HS_ENABLE_SHADER_WORKBENCH
