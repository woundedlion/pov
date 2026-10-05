/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file composed_effect.h
 * @brief Shared machinery for the composed-effect family: the parameter
 *        providers and present-only instance storage derived from a ranked
 *        stage Spec, with the engine's preset choreography and palette lifecycle.
 */

#include "animation/orientation.h"
#include "math/mobius.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <span>
#include <type_traits>

#include "color/effect_palette_recipes.h"
#include "control/choreography.h"
#include "color/palette_cycler.h"
#include "control/registry.h"
#include "engine/memory.h"
#include "render/scan.h"
#include "math/noise_field.h"
#include "render/pullback.h"
#include "render/pullback/runtime_seeds.h"
#include "render/pullback/composed_resources.h"

#if HS_ENABLE_TEST_HOOKS
namespace hs_test {
struct ComposedFrameWhiteBox;
}
#endif

namespace Pullback {

enum class SurfacePlacement : uint8_t { BEFORE_LENS, AFTER_LENS };

// The stage vocabulary a composed effect names, re-exported from the
// per-stage headers.
using Color::ColorParams;
using Color::HueMode;
using ValueCoverage::CutoutValueParams;
using ProjectionCoverage::EdgeValueParams;
using Lens::MobiusLensParams;
using Projection::ProjectionParams;
using Source::GridSourceParams;
using Source::LatticeSourceParams;
using Source::ProjectedNoiseSourceParams;
using Source::SphericalNoiseSourceParams;
using Source::SpiralSourceParams;
using Source::TwinWaveSourceParams;
using Surface::DirectSurfaceParams;
using Surface::PeriodicRippleParams;
using Surface::SurfaceNoiseParams;
using Transfer::IsoValueParams;
using Warp::AffineParams;
using Warp::MirrorParams;
using Warp::PolarParams;
using Warp::VectorNoiseParams;
using Warp::WaveShearParams;

/**
 * @brief The public frame context the pipeline's stages prepare from and
 *        shade against, resolved before the scan.
 * @details Read through the providers below; each stage's private prepared
 * state lives in the pipeline's per-frame instance instead. Its pointers
 * alias the runtime's persistent state and the palette cycler's current
 * bake, so a frame outlives only the draw_frame() call that built it.
 */
template <typename ParamsT> struct FrameState {
  /** Conjugate of the projection orientation; identity unless the effect sets
      `ANIMATED_PROJECTION`. */
  math::Quaternion projection_conjugate;
  /** Conjugate of the outer camera orientation. */
  math::Quaternion outer_conjugate;
  const BakedPalette *palette; /**< The cycler's current bake. */
  /** Hue-rotation LUT base; current only when hue_rotation_active(). */
  const Pixel *hue_rotation_lut;
  /** Hue-noise LUT base; current only under HueMode::NOISE with an active
      rotation. */
  const int8_t *hue_noise_lut;
  ParamsT params; /**< The frame's parameter values, already interpolated if a
                       preset transition is in flight. */
  /** Palette mapping weights, blended across a preset transition. */
  Pullback::Color::PaletteMappingWeights palette_mapping;
  ComposedDetail::ResourceStorage<ComposedDetail::ResourceFrame,
                                  typename ParamsT::ResourceTypes>
      resources;
  float palette_oscillation_phase; /**< Phase of the mapping wobble. */
};

/**
 * @brief Ties a pullback pipeline to one effect's frame state.
 * @details Composed effects render uninstrumented, so the binding pins
 * Pullback::NoInstrumentation.
 * @tparam FrameT The effect's FrameState specialization.
 */
template <typename FrameT> struct Binding {
  using FrameState = FrameT;
  using Instrumentation = Pullback::NoInstrumentation;
};

/** @brief Supplies the camera orientation to Pullback::Stage::Rotate. */
template <typename BindingT> struct OuterCameraProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  __attribute__((always_inline)) static const math::Quaternion &
  conjugate(const FrameState &frame) {
    return frame.outer_conjugate;
  }
};

/**
 * @brief Supplies the projection frame and its parameters to the
 *        Pullback::Projection policies.
 * @details Exposes the composed projection policies' accessors; an effect pays
 * only for the ones its chosen policy instantiates, so a projection that takes
 * no central meridian never reads that field.
 */
template <typename BindingT> struct ProjectionProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  __attribute__((always_inline)) static const math::Quaternion &
  conjugate(const FrameState &frame) {
    return frame.projection_conjugate;
  }
  __attribute__((always_inline)) static float
  singularity_fade(const FrameState &frame) {
    return frame.params.template get<"projection">().singularity_fade;
  }
  __attribute__((always_inline)) static float
  central_meridian(const FrameState &frame) {
    return frame.params.template get<"projection">().central_meridian;
  }
};

/** @brief Supplies the Mobius coefficients to Pullback::Lens::Mobius. */
template <typename BindingT, ResourceKey Key = "lens"> struct LensProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  __attribute__((always_inline)) static const math::MobiusParams &
  params(const FrameState &frame) {
    return frame.params.template get<Key>().mobius;
  }
};

/** @brief Supplies parameters, clock and noise of one planar warp instance. */
template <typename BindingT, ResourceKey Key, typename Family,
          bool TrackPath = false, ResourceKey SourceKey = "source">
struct WarpProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  __attribute__((always_inline)) static const auto &
  params(const FrameState &frame) {
    return frame.params.template get<Key>();
  }
  __attribute__((always_inline)) static auto prepare(const FrameState &frame) {
    using WarpT = std::remove_cvref_t<decltype(params(frame))>;
    if constexpr (std::is_same_v<WarpT, AffineParams>) {
      static_assert(
          requires {
            frame.params.template get<SourceKey>().lattice_cell_scale;
          }, "the affine warp stage translates in lattice cells and requires a "
             "LatticeSourceParams source");
      return Pullback::Warp::prepare(
          params(frame), phase(frame),
          frame.resources.template get<Key>().rotation,
          1.0f / frame.params.template get<SourceKey>().lattice_cell_scale);
    } else if constexpr (std::is_same_v<WarpT, PolarParams>) {
      return Pullback::NoPrepared{};
    } else {
      return Pullback::Warp::prepare(params(frame), phase(frame));
    }
  }
  __attribute__((always_inline)) static float phase(const FrameState &frame) {
    return frame.resources.template get<Key>().phase;
  }
  __attribute__((always_inline)) static const FastNoiseLite &
  noise(const FrameState &frame) {
    return *frame.resources.template get<Key>().noise;
  }
  __attribute__((always_inline)) static bool
  path_length_required(const FrameState &) {
    return TrackPath;
  }
};

/**
 * @brief Supplies the displacement field to the Pullback::Surface policies.
 * @tparam BindingT The effect's Binding.
 * @tparam TrackPath Whether the stage accumulates path length.
 * @pre The resource family is a displacement family.
 */
template <typename BindingT, typename Family, bool TrackPath = false,
          ResourceKey Key = "surface">
struct SurfaceProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  __attribute__((always_inline)) static const FastNoiseLite &
  noise(const FrameState &frame) {
    return *frame.resources.template get<Key>().noise;
  }
  __attribute__((always_inline)) static auto prepare(const FrameState &frame) {
    if constexpr (requires { frame.params.template get<Key>().direction; })
      return Pullback::Surface::prepare_direct(
          frame.resources.template get<Key>().phase,
          frame.params.template get<Key>().direction);
    else
      return Pullback::Surface::prepare(
          frame.resources.template get<Key>().phase);
  }
  __attribute__((always_inline)) static const auto &
  params(const FrameState &frame) {
    return frame.params.template get<Key>();
  }
  __attribute__((always_inline)) static float phase(const FrameState &frame) {
    using SurfaceParams = Family;
    if constexpr (std::is_same_v<SurfaceParams, PeriodicRippleParams>)
      return frame.resources.template get<Key>().phase /
             frame.params.template get<Key>().period;
    else
      return frame.resources.template get<Key>().phase;
  }
  __attribute__((always_inline)) static float scale(const FrameState &frame) {
    return frame.params.template get<Key>().scale;
  }
  __attribute__((always_inline)) static float
  strength(const FrameState &frame) {
    return frame.params.template get<Key>().strength;
  }
  __attribute__((always_inline)) static bool
  path_length_required(const FrameState &) {
    return TrackPath;
  }
};

/**
 * @brief Supplies the pattern and noise state to the Pullback::Source policies.
 * @details The pattern accessors read the
 * prepared phases, the noise accessors the NoiseSourceParams fields. Only the
 * accessors an effect's chosen source policy names are instantiated.
 */
template <typename BindingT, typename Family, ResourceKey Key = "source">
struct SourceProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  __attribute__((always_inline)) static const auto &
  params(const FrameState &frame) {
    return frame.params.template get<Key>();
  }
  __attribute__((always_inline)) static Pullback::Source::PreparedSource
  prepare(const FrameState &frame) {
    return Pullback::Source::prepare(
        frame.resources.template get<Key>().primary,
        frame.resources.template get<Key>().secondary,
        frame.resources.template get<Key>().angle);
  }
  __attribute__((always_inline)) static const FastNoiseLite &
  noise(const FrameState &frame) {
    return *frame.resources.template get<Key>().noise;
  }
  __attribute__((always_inline)) static float
  noise_scale(const FrameState &frame) {
    return frame.params.template get<Key>().noise_scale;
  }
  __attribute__((always_inline)) static float
  noise_time(const FrameState &frame) {
    return frame.resources.template get<Key>().noise_time;
  }
  __attribute__((always_inline)) static float
  noise_contrast(const FrameState &frame) {
    return frame.params.template get<Key>().noise_contrast;
  }
};

/**
 * @brief Supplies the value-family fields to the Pullback::Transfer and
 *        Pullback::ValueCoverage and ProjectionCoverage::EdgeFade policies.
 * @details Names the five shared value-family fields; an effect's material
 * stage instantiates only the accessors its transfer and coverage policies
 * call, so an IsoValueParams effect never touches `edge_width` and vice versa.
 */
template <typename BindingT, typename Family, ResourceKey Key = "value">
struct ValueProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  __attribute__((always_inline)) static float
  iso_level(const FrameState &frame) {
    return frame.params.template get<Key>().iso_level;
  }
  __attribute__((always_inline)) static float
  iso_width(const FrameState &frame) {
    return frame.params.template get<Key>().iso_width;
  }
  __attribute__((always_inline)) static float
  edge_width(const FrameState &frame) {
    return frame.params.template get<Key>().edge_width;
  }
  __attribute__((always_inline)) static float
  cutout_threshold(const FrameState &frame) {
    return frame.params.template get<Key>().cutout_threshold;
  }
  __attribute__((always_inline)) static float
  cutout_softness(const FrameState &frame) {
    return frame.params.template get<Key>().cutout_softness;
  }
};

/**
 * @brief Whether the colorizer samples the hue-rotation LUT this frame.
 * @details The runtime rebuilds the LUT on exactly this condition, so the two
 * sites cannot disagree about which frames leave it stale.
 */
template <HueMode HueV>
inline bool hue_rotation_active(const ColorParams &color) {
  return HueV != HueMode::NONE && color.hue_shift_amount != 0.0f;
}

/**
 * @brief Supplies the palette, mapping and hue state to
 *        Pullback::Color::GeneratedPalette.
 * @details Both LUT views carry their own active flag, so a stale LUT is never
 * sampled: the noise view additionally requires HueMode::NOISE.
 * @tparam BindingT The effect's Binding.
 * @tparam HueV Hue-rotation source reported to the color stage.
 * @tparam BrightnessV Brightness envelope reported to the color stage.
 */
template <typename BindingT, HueMode HueV,
          Pullback::Color::BrightnessEnvelope BrightnessV>
struct ColorProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  __attribute__((always_inline)) static Pullback::Color::PaletteMappingWeights
  mapping_weights(const FrameState &frame) {
    return frame.palette_mapping;
  }
  __attribute__((always_inline)) static float
  mapping_frequency(const FrameState &frame) {
    return frame.params.template get<"color">().mapping_frequency;
  }
  __attribute__((always_inline)) static float
  mapping_phase(const FrameState &frame) {
    return frame.params.template get<"color">().mapping_phase;
  }
  __attribute__((always_inline)) static float
  oscillation_depth(const FrameState &frame) {
    return frame.params.template get<"color">().phase_oscillation_depth;
  }
  __attribute__((always_inline)) static float
  oscillation_phase(const FrameState &frame) {
    return frame.palette_oscillation_phase;
  }
  __attribute__((always_inline)) static const BakedPalette &
  palette(const FrameState &frame) {
    return *frame.palette;
  }
  __attribute__((always_inline)) static Pullback::Color::HueMode
  hue_mode(const FrameState &) {
    return HueV;
  }
  __attribute__((always_inline)) static float
  hue_shift_amount(const FrameState &frame) {
    return frame.params.template get<"color">().hue_shift_amount;
  }
  __attribute__((always_inline)) static Pullback::Color::HueRotationLutView
  hue_rotation(const FrameState &frame) {
    return {frame.hue_rotation_lut,
            hue_rotation_active<HueV>(frame.params.template get<"color">())};
  }
  __attribute__((always_inline)) static Pullback::Color::HueNoiseLutView
  hue_noise(const FrameState &frame) {
    return {frame.hue_noise_lut, HueV == HueMode::NOISE &&
                                     hue_rotation_active<HueV>(
                                         frame.params.template get<"color">())};
  }
  __attribute__((always_inline)) static Pullback::Color::BrightnessEnvelope
  brightness_envelope(const FrameState &) {
    return BrightnessV;
  }
  __attribute__((always_inline)) static float
  brightness_bottom(const FrameState &frame) {
    return frame.params.template get<"color">().brightness_bottom;
  }
  __attribute__((always_inline)) static float
  brightness_top(const FrameState &frame) {
    return frame.params.template get<"color">().brightness_top;
  }
  __attribute__((always_inline)) static float
  opacity_low(const FrameState &frame) {
    return frame.params.template get<"color">().opacity_low;
  }
  __attribute__((always_inline)) static float
  opacity_high(const FrameState &frame) {
    return frame.params.template get<"color">().opacity_high;
  }
};

/**
 * @brief Interpolates one parameter family across a preset transition.
 * @details Driven by the family's field table; each field moves on the curve
 * its descriptor names, and every member the table does not cover, the Mobius
 * coefficients included, snaps to @p b once progress reaches 1.
 * @param a Value at progress 0.
 * @param b Value at progress 1.
 * @param t Progress fraction.
 * @return The interpolated family.
 */
template <Pullback::HasFields T>
inline T interpolate(const T &a, const T &b, float t) {
  return Pullback::Fields::interpolate(a, b, t);
}

inline MobiusLensParams interpolate(const MobiusLensParams &a,
                                    const MobiusLensParams &b, float t) {
  MobiusLensParams value;
  value.mobius = t < 1.0f ? a.mobius : b.mobius;
  return value;
}

/**
 * @brief Interpolates a whole parameter set, family by family.
 * @param from Parameters at progress 0.
 * @param to Parameters at progress 1.
 * @param progress Progress fraction, already eased by the caller.
 * @return The interpolated parameter set.
 */
template <typename... Resources>
inline ComposedDetail::ParameterSet<ComposedDetail::ResourceList<Resources...>>
interpolate(const ComposedDetail::ParameterSet<
                ComposedDetail::ResourceList<Resources...>> &from,
            const ComposedDetail::ParameterSet<
                ComposedDetail::ResourceList<Resources...>> &to,
            float progress) {
  return {ComposedDetail::ParameterBlock<Resources::KEY,
                                         typename Resources::Family>{
      interpolate(from.template get<Resources::KEY>(),
                  to.template get<Resources::KEY>(), progress)}...};
}

/**
 * @brief Whether every field of a parameter family is inside its authored
 *        range.
 * @details Driven by the family's field table, so the admissibility ranges a
 * restored snapshot must pass are the same descriptors the sliders register
 * with.
 */
template <Pullback::HasFields T> inline bool valid(const T &value) {
  return Pullback::Fields::valid(value);
}

inline bool valid(const MobiusLensParams &p) {
  const float values[] = {p.mobius.a.re, p.mobius.a.im, p.mobius.b.re,
                          p.mobius.b.im, p.mobius.c.re, p.mobius.c.im,
                          p.mobius.d.re, p.mobius.d.im};
  for (float value : values)
    if (!std::isfinite(value) ||
        fabsf(value) > MobiusLensParams::COEFFICIENT_LIMIT)
      return false;
  return MobiusLensParams::nondegenerate(p.mobius);
}

inline bool valid(const ColorParams &p) {
  return Pullback::Fields::valid(p) &&
         static_cast<uint8_t>(p.palette_mapping) <=
             static_cast<uint8_t>(Pullback::Color::PaletteMapping::REVERSE);
}

/**
 * @brief Whether every family of a parameter set is in range.
 * @return True only when every declared instance passes.
 */
template <typename... Resources>
inline bool valid(const ComposedDetail::ParameterSet<
                  ComposedDetail::ResourceList<Resources...>> &params) {
  bool result = true;
  params.visit([&]<typename Resource>(const auto &family) {
    result = valid(family) && result;
  });
  return result;
}

template <bool Enabled> struct OptionalNoise {};
template <> struct OptionalNoise<true> {
  FastNoiseLite noise;
};

/** @brief Hue-rotation LUT storage; empty when the effect never rotates hue. */
template <bool Enabled> struct OptionalHueRotationLut {};
template <> struct OptionalHueRotationLut<true> {
  std::array<Pixel, Pullback::Color::HueRotationLutView::SIZE> hue_rotation_lut;
  /** Palette bake the resident table was built from; 0 matches no bake,
      forcing the first build. */
  uint32_t hue_rotation_lut_bake = 0;
};

/** @brief Hue-noise LUT and the inputs it was baked from; empty unless the
    hue source is the noise field. */
template <bool Enabled> struct OptionalHueNoiseLut {};
template <> struct OptionalHueNoiseLut<true> {
  FastNoiseLite color_noise;
  std::array<int8_t, Pullback::Color::HueNoiseLutView::SIZE> hue_noise_lut;
  Pullback::Color::HueNoiseBakeCache hue_noise_bake;
};

/** @brief Projection-walk noise storage; empty when disabled. */
template <bool Enabled> struct ProjectionWalkNoise {};
template <> struct ProjectionWalkNoise<true> {
  FastNoiseLite projection_walk_noise;
};

/** @brief Persistent projection-walk state; empty when disabled. */
template <bool Enabled> struct ProjectionWalkState {
  math::Quaternion frame_conjugate() const { return math::Quaternion(); }
};
template <> struct ProjectionWalkState<true> {
  math::Orientation<> projection_walk;
  math::Quaternion projection_walk_previous;
  math::Quaternion projection_wander;
  math::Quaternion projection_conjugate;
  math::Quaternion base_orientation = Pullback::projection_base_orientation();
  float projection_spin = 0.0f;

  math::Quaternion frame_conjugate() const { return projection_conjugate; }
};

/** @brief Sphere-to-plane projection of a composed effect's Stage::Project. */
enum class ProjectionKind : uint8_t {
  STEREOGRAPHIC,
  GNOMONIC_FOLDED,
  EQUIRECTANGULAR,
  FOLDED_SINUSOIDAL
};

/** @brief Whether @p projection reads the central-meridian field. */
constexpr bool uses_central_meridian(ProjectionKind projection) {
  return projection == ProjectionKind::EQUIRECTANGULAR ||
         projection == ProjectionKind::FOLDED_SINUSOIDAL;
}

/** @brief Whether @p projection reads the singularity-fade field. Folded
    sinusoidal has no singular locus and returns fixed weights. */
constexpr bool uses_singularity_fade(ProjectionKind projection) {
  return projection != ProjectionKind::FOLDED_SINUSOIDAL;
}

/** @brief Optional transfer curve an effect's material stage composes. */
enum class TransferKind : uint8_t { NONE, ISO_CONTOUR };

/** @brief Optional value-dependent coverage stage after sampling. */
enum class FieldCoverageKind : uint8_t { NONE, VALUE_CUTOUT };

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

/**
 * @brief A complete composed effect: the shared lifecycle plus a pipeline
 *        declared by a ranked stage Spec.
 * @details The effect states its Spec and identity constants; parameter and
 * runtime storage derive from its pipeline providers. The shared
 * lifecycle — parameter registration, preset choreography, palette cycling,
 * camera walks and noise clocks — are assembled here. Required `Derived`
 * members are EFFECT_ID, PRESET_IDS, PARAMETER_SCHEMA_VERSION and
 * PRESET_DWELL_FRAMES. Presets and their departures resolve
 * through `preset(index)`, then `PRESETS`; only single-preset effects may fall
 * back to startup params.
 * `initial_params` is optional. Other optional members are `ANIMATED_MOBIUS`,
 * `CAMERA_SPIN_RATE` and an `after_composed_init()` hook; `WARP_NOISE_SEED` /
 * `SOURCE_NOISE_SEED` / `SURFACE_NOISE_SEED` are inherited members an effect
 * shadows to decorrelate its warp, source, or surface noise fields. A shade() shadow that forwards to
 * RenderPipeline::shade changes only the entry trampoline's placement; the
 * pipeline body remains in hot flash. Different body emission requires calling
 * RenderPipeline::evaluate(view, frame.ctx, frame.prepared) from the shadow.
 * Surface-noise effects conventionally wrap the sphere run (displacement, lens
 * and projection) in Stage::Placed<CodeEmission::OUT_OF_LINE_FLASH, ...>.
 *
 * `EFFECT_ID` is the registry identity; `PRESET_IDS` lists immutable preset
 * identities indexed by preset number; `PARAMETER_SCHEMA_VERSION` changes
 * with the Params layout to reject stale snapshots; `PRESET_DWELL_FRAMES`
 * gives the frames held before the next transition.
 * `DESCRIPTOR_DIGEST` and `PRESET_BANK_DIGEST` pin the pattern document's
 * canonical descriptor (excluding parameter units) and preset bank. The
 * product-group generator parity tests check them; the runtime never reads
 * them. The browser computes its own matching digests from pattern documents.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam Derived The effect class deriving from this base.
 * @tparam SpecT The effect's `Spec`.
 */
template <int W, int H, typename Derived, typename SpecT>
class ComposedEffect : public ChoreographedEffect<Derived, ParamsFor<SpecT>>,
                       private ProjectionWalkState<SpecT::ANIMATED_PROJECTION> {
  using ParamsT = ParamsFor<SpecT>;
  using Choreography = ChoreographedEffect<Derived, ParamsT>;
  friend Choreography;
  static constexpr PaletteHarmony Harmony = SpecT::HARMONY;
  static constexpr HueMode HueV = SpecT::HUE;
  static constexpr Color::BrightnessEnvelope BrightnessV = SpecT::BRIGHTNESS;
  static constexpr bool AnimatedProjection = SpecT::ANIMATED_PROJECTION;
#if HS_ENABLE_TEST_HOOKS
  friend struct hs_test::ComposedFrameWhiteBox;
#endif

public:
  using Params = ParamsT;
  using Spec = SpecT;
  static_assert(
      Spec::TRANSFER != TransferKind::ISO_CONTOUR || requires(Params p) {
        p.template get<"value">().iso_level;
        p.template get<"value">().iso_width;
      }, "iso-contour transfer requires iso-level and iso-width");
  static_assert(
      Spec::COVERAGE != ProjectionCoverageMode::EDGE_FADE ||
          requires(Params p) { p.template get<"value">().edge_width; },
      "edge-fade coverage requires edge-width");
  static_assert(
      Spec::FIELD_COVERAGE != FieldCoverageKind::VALUE_CUTOUT ||
          requires(Params p) {
            p.template get<"value">().cutout_threshold;
            p.template get<"value">().cutout_softness;
          },
      "value-cutout coverage requires threshold and softness");
  using FrameState = Pullback::FrameState<ParamsT>;
  using Binding = Pullback::Binding<FrameState>;
  static constexpr bool ANIMATED_PROJECTION = AnimatedProjection;
  template <ResourceKind Kind> static consteval bool has_noise() {
    bool result = false;
    Params{}.visit([&]<typename Resource>(const auto &) {
      if constexpr (Resource::KIND == Kind &&
                    ComposedDetail::RESOURCE_NOISE<Resource>)
        result = true;
    });
    return result;
  }
  static constexpr bool HAS_SURFACE_NOISE = has_noise<ResourceKind::SURFACE>();
  static constexpr bool HAS_WARP_NOISE = has_noise<ResourceKind::WARP>();
  static constexpr bool HAS_SOURCE_NOISE = has_noise<ResourceKind::SOURCE>();
  using RenderPipeline = typename SpecT::template Pipeline<Binding>;
  using Metadata =
      ComposedDetail::PipelineMetadata<Spec, Binding, RenderPipeline>;
  static_assert(Metadata::PATH_TRACKED == (Spec::HUE == HueMode::PATH_LENGTH),
                "path length hue metadata must match pipeline tracking");
  static_assert(Metadata::COVERAGE_MATCHES,
                "projection coverage metadata must match the pipeline");
  static_assert(Metadata::LENS_MATCHES,
                "lens policy metadata must match the pipeline");
  static_assert(Metadata::PROJECTION_MATCHES,
                "projection metadata must match the pipeline");
  static_assert(Metadata::SURFACE_PLACEMENT_MATCHES,
                "surface placement metadata must match the pipeline");
  static_assert(
      RenderPipeline::template any_stage<ComposedDetail::IsTransferStage> ==
          (Spec::TRANSFER != TransferKind::NONE),
      "transfer metadata must match the pipeline");
  static_assert(
      RenderPipeline::template any_stage<ComposedDetail::IsCoverageStage> ==
          (Spec::FIELD_COVERAGE != FieldCoverageKind::NONE),
      "field coverage metadata must match the pipeline");
  using Frame = typename RenderPipeline::Frame;

  /** @brief Per-field noise seeds; an effect shadows one to decorrelate its
      field from the shared spatial phase. */
  static constexpr int32_t WARP_NOISE_SEED = EFFECT_NOISE_SEED;
  static constexpr int32_t SOURCE_NOISE_SEED = EFFECT_NOISE_SEED;
  static constexpr int32_t SURFACE_NOISE_SEED = EFFECT_NOISE_SEED;

  /** @brief Constructs the effect at W x H with the POV column strobe on. */
  HS_COLD_MEMBER ComposedEffect() : Choreography(W, H, {.strobe = true}) {}

  /**
   * @brief Claims the runtime's persistent storage and registers the effect.
   * @details Every allocation this makes comes from `persistent_arena`, so the
   * effect heap-allocates nothing after init. `Derived::after_composed_init()`
   * runs last, once the parameters, palette cycler and camera walks are all
   * live, which is what lets it adjust the dwell or start a timeline
   * animation.
   */
  HS_COLD_MEMBER void init() override {
    this->begin_choreography();
    state = persistent_arena.make<State>();
    use_parameter_storage(persistent_arena,
                          persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);
    if constexpr (HueV == HueMode::NOISE)
      Pullback::init_effect_noise(state->color_noise, HUE_NOISE_SEED);
    params.visit([&]<typename Resource>(auto &) {
      if constexpr (ComposedDetail::RESOURCE_NOISE<Resource>) {
        constexpr int32_t SEED = Resource::KIND == ResourceKind::WARP
                                     ? Derived::WARP_NOISE_SEED
                                 : Resource::KIND == ResourceKind::SOURCE
                                     ? Derived::SOURCE_NOISE_SEED
                                     : Derived::SURFACE_NOISE_SEED;
        constexpr int32_t INSTANCE_SEED = [] {
          if constexpr (ComposedDetail::qualified<
                            Resource, typename Params::ResourceTypes>()) {
            constexpr uint32_t HASH = fnv1a(Resource::KEY.view());
            return static_cast<int32_t>(static_cast<uint32_t>(SEED) ^ HASH);
          } else {
            return SEED;
          }
        }();
        Pullback::init_effect_noise(
            state->resources.template get<Resource::KEY>().noise,
            INSTANCE_SEED);
      }
    });
    palette_cycler.init_generated(persistent_arena, next_palette, this, 0, 600,
                                  math::ease_in_out_sin);
    if constexpr (AnimatedProjection)
      timeline.add(0, Animation::RandomWalk<W>(
                          this->projection_walk, math::UP,
                          state->projection_walk_noise,
                          typename Animation::RandomWalk<W>::Options{},
                          PROJECTION_WALK_SEED));
    timeline.add(0, Animation::RandomWalk<W>(
                        outer_walk, math::UP, state->outer_walk_noise,
                        typename Animation::RandomWalk<W>::Options{},
                        CAMERA_WALK_SEED));
    register_parameters();
    if constexpr (requires(Derived &effect) { effect.after_composed_init(); })
      static_cast<Derived &>(*this).after_composed_init();
  }

  /**
   * @brief Advances the frame clocks, resolves the frame state and shades.
   * @details The order is load-bearing: the runtime's clocks, camera walks and
   * palette all step before prepare_frame() snapshots them, so the whole scan
   * shades from one consistent frame. The preset-transition lerp steps with
   * the timeline, so its writes land before the snapshot too.
   */
  HS_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(fx_timeline_step);
      timeline.step(canvas);
    }
    {
      HS_PROFILE(fx_advance);
      step_choreography();
      advance_runtime();
      update_spatial_frames();
      update_palette_chroma();
      palette_cycler.step();
    }
    const typename Derived::RenderPipeline::Frame frame =
        Derived::RenderPipeline::prepare(prepare_frame());
    {
      HS_PROFILE(fx_shader_draw);
      Scan::Shader::draw_cached<W, H, 1>(canvas,
                                         [&frame](const math::Vector &view) {
                                           return Derived::shade(view, frame);
                                         });
    }
  }

  /**
   * @brief Shades one pixel; forwards to RenderPipeline::shade.
   * @param view Unit view direction for the pixel.
   * @param frame Per-frame transforms, params and LUTs from the runtime.
   */
  static HS_O3_FN Color4 shade(const math::Vector &view, const Frame &frame) {
    return RenderPipeline::shade(view, frame);
  }

  /** @brief Whether a parameter set is admissible, family by family. */
  static bool valid_params(const Params &params) {
    return Pullback::valid(params);
  }
#if HS_ENABLE_PARAM_GUI_BRIDGE
  const char *parameter_warning(const char *name) const override {
    return refused_name != nullptr && std::strcmp(name, refused_name) == 0
               ? Lens::MobiusLensParams::DEGENERATE_WARNING
               : nullptr;
  }
#endif
#if HS_ENABLE_TEST_HOOKS
  FrameState frame_for_test() { return prepare_frame(); }
#endif

protected:
  /** Descriptors the arena-backed parameter array holds; every slider an effect
      registers has to fit. */
  static constexpr size_t PARAM_CAPACITY =
      ComposedDetail::parameter_count<SpecT>(typename Params::ResourceTypes{});

  using Choreography::anims_paused;
  using Choreography::params;

#if HS_ENABLE_PARAM_GUI_BRIDGE
  bool parameter_write_admitted(const ParamDef &parameter,
                                float value) override {
    bool admitted = true;
    params.visit([&]<typename Resource>(const auto &family) {
      if constexpr (Resource::KIND == ResourceKind::LENS) {
        const uintptr_t begin = reinterpret_cast<uintptr_t>(&family);
        const uintptr_t target = reinterpret_cast<uintptr_t>(parameter.target);
        if (target >= begin && target - begin < sizeof(family)) {
          auto candidate = family;
          ParamDef proposed = parameter;
          proposed.target =
              reinterpret_cast<unsigned char *>(&candidate) + (target - begin);
          this->write_parameter_unchecked(proposed, value);
          admitted = Pullback::valid(candidate);
        }
      }
    });
    if (!admitted) {
      refused_name = parameter.name;
      return false;
    }
    refused_name = nullptr;
    return true;
  }

  const char *refused_name = nullptr;
#endif

  using Choreography::step_choreography;
  using Choreography::register_animated_param;
  using Choreography::timeline;
  using Choreography::transition;
  using Choreography::use_parameter_storage;

  template <typename Resource>
  HS_COLD_MEMBER const char *resource_parameter_name(const char *name) {
    constexpr bool QUALIFY =
        ComposedDetail::qualified<Resource, typename Params::ResourceTypes>();
    if constexpr (!QUALIFY)
      return name;
    else {
      constexpr auto KEY = Resource::KEY.view();
      const size_t length = std::strlen(name);
      char *qualified =
          persistent_arena.allocate_n<char>(KEY.size() + length + 2);
      std::memcpy(qualified, KEY.data(), KEY.size());
      qualified[KEY.size()] = '.';
      std::memcpy(qualified + KEY.size() + 1, name, length + 1);
      return qualified;
    }
  }

  /** @brief Registers the named descriptors of one independent instance. */
  template <typename Resource, typename T>
  HS_COLD_MEMBER void register_fields(T &family) {
    for (const auto &field : T::FIELDS)
      if (field.name != nullptr && field_gate_open(field.gate))
        this->register_param(resource_parameter_name<Resource>(field.name),
                             &(family.*field.member), field.description().spec);
  }

  /** @brief Adopts a snap target and re-derives the palette mapping weights. */
  HS_COLD_MEMBER void adopt_params(const Params &target) {
    params = target;
    palette_mapping = Pullback::Color::PaletteMappingWeights::single(
        target.template get<"color">().palette_mapping);
  }

  /** @brief Adopts an automatic target while retaining animated lens coefficients. */
  HS_COLD_MEMBER void finish_blend(const Params &target)
    requires(requires { Derived::ANIMATED_MOBIUS; } && Derived::ANIMATED_MOBIUS)
  {
    params.visit([&]<typename Resource>(auto &family) {
      if constexpr (Resource::KIND == ResourceKind::LENS) {
        const auto MOBIUS = family.mobius;
        family = target.template get<Resource::KEY>();
        family.mobius = MOBIUS;
      } else {
        family = target.template get<Resource::KEY>();
      }
    });
    palette_mapping = Pullback::Color::PaletteMappingWeights::single(
        target.template get<"color">().palette_mapping);
  }

  HS_COLD_MEMBER void parameter_written() override {
    Choreography::parameter_written();
    palette_mapping = Pullback::Color::PaletteMappingWeights::single(
        params.template get<"color">().palette_mapping);
  }

  /** @brief Captures the palette-mapping endpoints of an arming crossfade. */
  HS_COLD_MEMBER void transition_armed(const Params &target) {
    mapping_from = palette_mapping;
    mapping_to = Pullback::Color::PaletteMappingWeights::single(
        target.template get<"color">().palette_mapping);
  }

  /**
   * @brief Writes the interpolated parameters of an in-flight transition.
   * @details Under `Derived::ANIMATED_MOBIUS` the live lens coefficients are
   * carried across the interpolation, so the transition does not overwrite the
   * timeline animation driving them.
   */
  HS_COLD_MEMBER void blend_params(float progress) {
    const auto before = params;
    params = Pullback::interpolate(transition.from, transition.to, progress);
    if constexpr (requires { Derived::ANIMATED_MOBIUS; })
      if constexpr (Derived::ANIMATED_MOBIUS)
        params.visit([&]<typename Resource>(auto &family) {
          if constexpr (Resource::KIND == ResourceKind::LENS)
            family.mobius = before.template get<Resource::KEY>().mobius;
        });
    palette_mapping = Pullback::Color::PaletteMappingWeights::lerp(
        mapping_from, mapping_to, progress);
  }

  /**
   * @brief Starts a repeating circular Mobius warp on the lens parameters.
   * @details Pausable, so the warp freezes with the rest of the parameter
   * animations. Compiles to nothing for an effect whose lens family carries no
   * Mobius coefficients.
   * @param scale Radius of the circular warp.
   * @param duration Frames per revolution.
   */
  template <ResourceKey Key = "lens">
  HS_COLD_MEMBER void start_mobius_animation(float scale, int duration) {
    if constexpr (Params::template HAS<Key>)
      timeline.add_pausable(
          0,
          Animation::MobiusWarpCircular(params.template get<Key>().mobius,
                                        scale, duration, true),
          &anims_paused);
  }

private:
  struct State : ProjectionWalkNoise<AnimatedProjection>,
                 OptionalHueRotationLut<HueV != HueMode::NONE>,
                 OptionalHueNoiseLut<HueV == HueMode::NOISE> {
    ComposedDetail::ResourceStorage<ComposedDetail::ResourceNoise,
                                    typename Params::ResourceTypes>
        resources;
    FastNoiseLite outer_walk_noise;
  };

  // init() takes the palette cycler's generated arena, the external ParamDef
  // array, one State and qualified parameter names from the persistent arena.
  static constexpr size_t FOOTPRINT_BYTES =
      PaletteCycler::generated_arena_bytes() +
      PARAM_CAPACITY * sizeof(ParamDef) + sizeof(State) + alignof(State) +
      alignof(ParamDef) +
      ComposedDetail::parameter_name_bytes<SpecT>(
          typename Params::ResourceTypes{});
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
      "Pullback::ComposedEffect persistent footprint exceeds the default "
      "partition");

  /** @brief Whether a gated field's slider exists for this effect. */
  static constexpr bool field_gate_open(Pullback::FieldGate gate) {
    return ComposedDetail::field_gate_open<SpecT>(gate);
  }

  template <typename Resource, typename T>
  HS_COLD_MEMBER void register_warp_fields(T &warp) {
    static_assert(T::FIELDS[0].member == &T::speed &&
                      T::FIELDS[0].name == nullptr,
                  "warp speed must be the first unnamed descriptor");
    constexpr const char *SPEED_NAME =
        ComposedDetail::warp_speed_name(Resource::KEY.view()).data();
    this->register_param(resource_parameter_name<Resource>(SPEED_NAME),
                         &warp.speed, T::FIELDS[0].description().spec);
    register_fields<Resource>(warp);
  }

  template <typename Resource, typename T>
  HS_COLD_MEMBER void register_lens_fields(T &lens) {
    constexpr math::Complex math::MobiusParams::*COEFFICIENTS[] = {
        &math::MobiusParams::a, &math::MobiusParams::b, &math::MobiusParams::c,
        &math::MobiusParams::d};
    constexpr float math::Complex::*CHANNELS[] = {&math::Complex::re,
                                                  &math::Complex::im};
    constexpr float LIMIT = MobiusLensParams::COEFFICIENT_LIMIT;
    for (size_t i = 0; i < std::size(COEFFICIENTS); ++i)
      for (size_t j = 0; j < std::size(CHANNELS); ++j)
        register_animated_param(
            resource_parameter_name<Resource>(
                ComposedDetail::MOBIUS_PARAM_NAMES[i * 2 + j]),
            &((lens.mobius.*COEFFICIENTS[i]).*CHANNELS[j]), -LIMIT, LIMIT);
  }

  template <float Color::ColorControls::*Member>
  static consteval const Field<ColorParams> &color_descriptor() {
    constexpr size_t index = [] {
      for (size_t i = 0; i < ColorParams::FIELDS.size(); ++i)
        if (ColorParams::FIELDS[i].member == Member)
          return i;
      return ColorParams::FIELDS.size();
    }();
    static_assert(index < ColorParams::FIELDS.size(),
                  "ColorParams member has no field descriptor");
    return ColorParams::FIELDS[index];
  }

  template <float Color::ColorControls::*Member>
  HS_COLD_MEMBER void register_color_field(const char *name) {
    constexpr const auto &field = color_descriptor<Member>();
    register_animated_param(name, &(params.template get<"color">().*Member),
                            field.min, field.max);
  }

  template <ResourceKind Kind> HS_COLD_MEMBER void register_resource_kind() {
    params.visit([&]<typename Resource>(auto &family) {
      if constexpr (Resource::KIND == Kind) {
        if constexpr (Kind == ResourceKind::WARP)
          register_warp_fields<Resource>(family);
        else if constexpr (Kind == ResourceKind::LENS)
          register_lens_fields<Resource>(family);
        else
          register_fields<Resource>(family);
      }
    });
  }

  HS_COLD_MEMBER void register_parameters() {
    register_resource_kind<ResourceKind::SOURCE>();
    register_resource_kind<ResourceKind::PROJECTION>();
    register_resource_kind<ResourceKind::SURFACE>();
    register_resource_kind<ResourceKind::WARP>();
    register_resource_kind<ResourceKind::VALUE>();
    register_resource_kind<ResourceKind::LENS>();
    register_color_field<&ColorParams::palette_chroma>("Palette Chroma");
    register_animated_param(
        "Palette Mapping", &params.template get<"color">().palette_mapping,
        PALETTE_MAPPING_OPTIONS, PALETTE_MAPPING_EXPORT_OPTIONS,
        std::size(PALETTE_MAPPING_OPTIONS));
    register_color_field<&ColorParams::mapping_frequency>("Mapping Frequency");
    register_color_field<&ColorParams::mapping_phase>("Mapping Phase");
    register_color_field<&ColorParams::phase_oscillation_depth>(
        "Phase Oscillation Depth");
    register_color_field<&ColorParams::phase_oscillation_speed>(
        "Phase Oscillation Speed");
    if constexpr (BrightnessV != Pullback::Color::BrightnessEnvelope::NONE) {
      register_color_field<&ColorParams::brightness_bottom>(
          "Brightness Bottom");
      register_color_field<&ColorParams::brightness_top>("Brightness Top");
    }
    register_color_field<&ColorParams::opacity_low>("Opacity at Value 0");
    register_color_field<&ColorParams::opacity_high>("Opacity at Value 1");
    if constexpr (HueV != HueMode::NONE)
      register_color_field<&ColorParams::hue_shift_amount>("Hue Shift Amount");
    if constexpr (HueV == HueMode::NOISE) {
      register_color_field<&ColorParams::hue_noise_scale>("Hue Noise Scale");
      register_color_field<&ColorParams::hue_noise_speed>("Hue Noise Speed");
    }
  }

  /**
   * @brief Steps every phase clock the effect's parameter families define.
   * @details Each clock is compiled in only when its field exists, so an effect
   * pays for the clocks of its declared instances. Affine rotations are
   * accumulated independently for each affine warp.
   */
  HS_COLD_MEMBER void advance_runtime() {
    params.visit([&]<typename Resource>(const auto &family) {
      auto &clock = clocks.template get<Resource::KEY>();
      if constexpr (Resource::KIND == ResourceKind::SOURCE) {
        Source::advance_clocks(family, clock.primary, clock.secondary,
                               clock.angle);
        if constexpr (requires { family.noise_time_rate; })
          clock.noise_time =
              math::wrap_t(clock.noise_time + family.noise_time_rate);
      } else if constexpr (Resource::KIND == ResourceKind::SURFACE) {
        if constexpr (std::is_same_v<typename Resource::Family,
                                     PeriodicRippleParams>)
          clock.phase = fmodf(clock.phase + 1.0f, family.period);
        else
          clock.phase = math::wrap_t(clock.phase + family.speed);
      } else if constexpr (Resource::KIND == ResourceKind::WARP) {
        if constexpr (std::is_same_v<typename Resource::Family, AffineParams>)
          clock.rotation = math::TWO_PI_F *
                           math::wrap_t((clock.rotation +
                                         family.speed * family.rotation_rate) /
                                        math::TWO_PI_F);
        clock.phase = math::wrap_t(clock.phase + family.speed);
      }
    });
    if constexpr (AnimatedProjection)
      this->projection_spin = fmodf(
          this->projection_spin + params.template get<"projection">().spin_rate,
          math::TWO_PI_F);
    if constexpr (requires { Derived::CAMERA_SPIN_RATE; })
      camera_spin =
          fmodf(camera_spin + Derived::CAMERA_SPIN_RATE, math::TWO_PI_F);
    if constexpr (HueV == HueMode::NOISE)
      hue_noise_phase = math::wrap_t(
          hue_noise_phase + params.template get<"color">().hue_noise_speed);
    palette_oscillation_phase =
        math::wrap_t(palette_oscillation_phase +
                     params.template get<"color">().phase_oscillation_speed);
  }

  // Rotation samples are eased within each walk step; chain walks apply the
  // recurrence directly and use different seeds. Nonzero wander differs.
  HS_COLD_MEMBER void update_spatial_frames() {
    // prepare_frame() reads projection_conjugate only for an animated
    // projection.
    if constexpr (AnimatedProjection) {
      const math::Quaternion projection = this->projection_walk.get();
      const math::Quaternion projection_delta =
          projection * this->projection_walk_previous.conjugate();
      this->projection_walk_previous = projection;
      this->projection_wander =
          (math::scaled_rotation_delta(
               projection_delta.normalized(),
               params.template get<"projection">().wander) *
           this->projection_wander)
              .normalized();
      this->projection_conjugate =
          (math::make_rotation(math::Y_AXIS, this->projection_spin) *
           this->base_orientation * this->projection_wander)
              .conjugate();
    }
    const math::Quaternion outer = outer_walk.get();
    const math::Quaternion outer_delta =
        outer * outer_walk_previous.conjugate();
    outer_walk_previous = outer;
    outer_wander = (math::scaled_rotation_delta(
                        outer_delta.normalized(),
                        params.template get<"projection">().camera_wander) *
                    outer_wander)
                       .normalized();
    if constexpr (requires { Derived::CAMERA_SPIN_RATE; })
      outer_conjugate =
          (math::make_rotation(math::Y_AXIS, camera_spin) * outer_wander)
              .conjugate();
    else
      outer_conjugate = outer_wander.conjugate();
  }

  /** @brief The hue-rotation LUT base, or null when no hue mode reads one. */
  const Pixel *hue_rotation_lut_data() const {
    if constexpr (HueV != HueMode::NONE)
      return state->hue_rotation_lut.data();
    else
      return nullptr;
  }

  /** @brief The hue-noise LUT base, or null outside HueMode::NOISE. */
  const int8_t *hue_noise_lut_data() const {
    if constexpr (HueV == HueMode::NOISE)
      return state->hue_noise_lut.data();
    else
      return nullptr;
  }

  /**
   * @brief Bakes the frame's LUTs and snapshots everything the scan reads.
   * @details Both LUT builds are gated on hue_rotation_active(), the flag the
   * returned frame hands the color stage; under it the hue-noise LUT is rebuilt
   * only when its scale or phase moved, and the hue-rotation LUT only when the
   * palette cycler rebaked its display LUT.
   * @return The frame state for this draw, valid until the next draw_frame().
   */
  HS_COLD_MEMBER FrameState prepare_frame() {
    HS_PROFILE(fx_prepare_frame);
    if constexpr (HueV == HueMode::NOISE) {
      if (hue_rotation_active<HueV>(params.template get<"color">())) {
        state->hue_noise_bake.refresh(
            state->hue_noise_lut, state->color_noise,
            params.template get<"color">().hue_noise_scale, hue_noise_phase);
      }
    }
    if constexpr (HueV != HueMode::NONE)
      if (hue_rotation_active<HueV>(params.template get<"color">()) &&
          state->hue_rotation_lut_bake != palette_cycler.bake_generation()) {
        Pullback::Color::prepare_hue_rotation_lut(
            std::span<Pixel, Pullback::Color::HueRotationLutView::SIZE>(
                state->hue_rotation_lut),
            palette_cycler.palette());
        state->hue_rotation_lut_bake = palette_cycler.bake_generation();
      }
    FrameState frame{.projection_conjugate = this->frame_conjugate(),
                     .outer_conjugate = outer_conjugate,
                     .palette = &palette_cycler.palette(),
                     .hue_rotation_lut = hue_rotation_lut_data(),
                     .hue_noise_lut = hue_noise_lut_data(),
                     .params = params,
                     .palette_mapping = palette_mapping,
                     .resources = {},
                     .palette_oscillation_phase = palette_oscillation_phase};
    params.visit([&]<typename Resource>(const auto &) {
      auto &resource = frame.resources.template get<Resource::KEY>();
      static_cast<ComposedDetail::ResourceClock<Resource> &>(resource) =
          clocks.template get<Resource::KEY>();
      if constexpr (ComposedDetail::RESOURCE_NOISE<Resource>)
        resource.noise = &state->resources.template get<Resource::KEY>().noise;
      if constexpr (Resource::KIND == ResourceKind::LENS)
        if (!Pullback::valid(frame.params.template get<Resource::KEY>()))
          frame.params.template get<Resource::KEY>().mobius = {};
    });
    return frame;
  }

  HS_COLD_MEMBER void update_palette_chroma() {
    if (palette_chroma == params.template get<"color">().palette_chroma)
      return;
    palette_chroma = params.template get<"color">().palette_chroma;
    palette_cycler.set_generated_chroma(palette_chroma);
  }

  /**
   * @brief Generator the palette cycler calls for each palette in the cycle.
   * @details Advances the hue on every palette after the first, so the cycle
   * opens on `palette_hue`, which starts at 0 and moves only here.
   * @param context The ComposedEffect instance, as registered with the cycler.
   * @param sequence Zero-based index of the palette being generated.
   * @param out Receives the recipe to bake.
   */
  static void next_palette(void *context, uint32_t sequence,
                           GenerativePalette &out) {
    ComposedEffect &effect = *static_cast<ComposedEffect *>(context);
    if (sequence > 0)
      effect.palette_hue += 159;
    out = GenerativePalette{PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, Harmony, AxisCurve::ASCENDING,
        PaletteRecipes::hue_turns(effect.palette_hue),
        effect.params.template get<"color">().palette_chroma)};
  }

  static constexpr const char *PALETTE_MAPPING_OPTIONS[] = {
      "Cup", "Bell", "Linear", "Reverse"};
  static constexpr const char *PALETTE_MAPPING_EXPORT_OPTIONS[] = {
      "Pullback::Color::PaletteMapping::CUP",
      "Pullback::Color::PaletteMapping::BELL",
      "Pullback::Color::PaletteMapping::LINEAR",
      "Pullback::Color::PaletteMapping::REVERSE"};

  State *state = nullptr;
  Pullback::Color::PaletteMappingWeights palette_mapping =
      Pullback::Color::PaletteMappingWeights::single(
          params.template get<"color">().palette_mapping);
  Pullback::Color::PaletteMappingWeights mapping_from;
  Pullback::Color::PaletteMappingWeights mapping_to;
  math::Orientation<> outer_walk;
  math::Quaternion outer_walk_previous;
  math::Quaternion outer_wander;
  math::Quaternion outer_conjugate;
  ComposedDetail::ResourceStorage<ComposedDetail::ResourceClock,
                                  typename Params::ResourceTypes>
      clocks;
  float camera_spin = 0.0f;
  float hue_noise_phase = 0.0f;
  float palette_oscillation_phase = 0.0f;
  float palette_chroma = -1.0f;
  uint32_t palette_hue = 0;
  PaletteCycler palette_cycler;
};

} // namespace Pullback
