/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file composed_effect.h
 * @brief Shared machinery for the composed-effect family: the parameter
 *        families an effect assembles into a `Params`, the frame-state
 *        providers its pullback pipeline reads, and the `ComposedEffect` base
 *        that assembles the pipeline over the engine's preset choreography.
 */

#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <type_traits>

#include "color/effect_palette_recipes.h"
#include "engine/control/choreography.h"
#include "engine/engine.h"
#include "math/noise_field.h"
#include "render/pullback.h"

namespace Pullback {

// The stage vocabulary a composed effect names, re-exported from the
// per-stage headers.
using Color::ColorParams;
using Color::HueMode;
using Coverage::EdgeValueParams;
using Lens::MobiusLensParams;
using Lens::NoLensParams;
using Projection::ProjectionParams;
using Source::GridSourceParams;
using Source::LatticeSourceParams;
using Source::NoiseSourceParams;
using Source::SpiralSourceParams;
using Source::TwinWaveSourceParams;
using Surface::DirectSurfaceParams;
using Surface::NoSurfaceParams;
using Surface::SurfaceNoiseParams;
using Transfer::IsoValueParams;
using Transfer::LinearValueParams;
using Warp::AffineParams;
using Warp::MirrorParams;
using Warp::NoWarpParams;
using Warp::PolarParams;
using Warp::VectorNoiseParams;
using Warp::WaveShearParams;

/**
 * @brief One composed effect's complete parameter set, one family per
 *        pipeline stage.
 * @details The member typedefs are the runtime's dispatch keys: parameter
 * registration, per-frame preparation and the warp/surface clocks all select
 * behavior by comparing them against the concrete families.
 * @tparam SourceT Source family, e.g. GridSourceParams.
 * @tparam OuterWarpT Family for the first planar warp slot.
 * @tparam InnerWarpT Family for the second planar warp slot.
 * @tparam LensT Lens family.
 * @tparam ValueT Value family read by the material stage.
 * @tparam SurfaceT Surface-displacement family.
 */
template <typename SourceT, typename OuterWarpT, typename InnerWarpT,
          typename LensT = NoLensParams, typename ValueT = LinearValueParams,
          typename SurfaceT = NoSurfaceParams>
struct Params {
  using source_type = SourceT;
  using outer_warp_type = OuterWarpT;
  using inner_warp_type = InnerWarpT;
  using lens_type = LensT;
  using value_type = ValueT;
  using surface_type = SurfaceT;
  SourceT source;
  ProjectionParams projection;
  OuterWarpT outer_warp;
  InnerWarpT inner_warp;
  SurfaceT surface;
  LensT lens;
  ValueT value;
  ColorParams color;
};

/**
 * @brief The public frame context the pipeline's stages prepare from and
 *        shade against, resolved before the scan.
 * @details Read through the providers below; each stage's private prepared
 * state lives in the pipeline's per-frame instance instead. Its pointers
 * alias the runtime's persistent state and the palette cycler's current
 * bake, so a frame outlives only the draw_frame() call that built it.
 * @tparam ParamsT The effect's `Pullback::Params` specialization.
 */
template <typename ParamsT> struct FrameState {
  /** Conjugate of the projection orientation; identity unless the effect sets
      `ANIMATED_PROJECTION`. */
  Quaternion projection_conjugate;
  /** Conjugate of the outer camera orientation. */
  Quaternion outer_conjugate;
  const FastNoiseLite *outer_noise;   /**< Null unless `HasOuterNoise`. */
  const FastNoiseLite *source_noise;  /**< Null unless `HasSourceNoise`. */
  const FastNoiseLite *surface_noise; /**< Null unless the effect displaces. */
  const BakedPalette *palette;        /**< The cycler's current bake. */
  /** Hue-rotation LUT base; current only when hue_rotation_active(). */
  const Pixel *hue_rotation_lut;
  /** Hue-noise LUT base; current only under HueMode::NOISE with an active
      rotation. */
  const int8_t *hue_noise_lut;
  ParamsT params; /**< The frame's parameter values, already interpolated if a
                       preset transition is in flight. */
  /** Palette mapping weights, blended across a preset transition. */
  Pullback::Color::PaletteMappingWeights palette_mapping;
  float source_primary;            /**< Source primary phase. */
  float source_secondary;          /**< Source secondary phase. */
  float source_angle;              /**< Source rotation, in radians. */
  float outer_phase;               /**< First warp slot's phase. */
  float inner_phase;               /**< Second warp slot's phase. */
  float outer_rotation;            /**< First slot's accumulated rotation. */
  float source_noise_time;         /**< Source noise time coordinate. */
  float surface_phase;             /**< Displacement-loop phase. */
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

/** @brief Supplies the camera orientation to Pullback::Stage::OuterCamera. */
template <typename BindingT> struct OuterCameraProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static const Quaternion &conjugate(const FrameState &frame) {
    return frame.outer_conjugate;
  }
};

/**
 * @brief Supplies the projection frame and its parameters to the
 *        Pullback::Projection policies.
 * @details Exposes every accessor the projection policies name; an effect pays
 * only for the ones its chosen policy instantiates, so a projection that takes
 * no central meridian never reads that field.
 */
template <typename BindingT> struct ProjectionProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static const Quaternion &conjugate(const FrameState &frame) {
    return frame.projection_conjugate;
  }
  static float pole_fade(const FrameState &frame) {
    return frame.params.projection.pole_fade;
  }
  static float central_meridian(const FrameState &frame) {
    return frame.params.projection.central_meridian;
  }
};

/** @brief Supplies the Mobius coefficients to Pullback::Lens::Mobius. */
template <typename BindingT> struct LensProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static const MobiusParams &params(const FrameState &frame) {
    return frame.params.lens.mobius;
  }
};

/**
 * @brief Supplies one planar warp slot to the Pullback::Warp policies.
 * @details `noise()` returns the outer noise field for either slot, so a
 * noise-driven warp in the inner slot still requires the runtime's
 * `HasOuterNoise`.
 * @tparam BindingT The effect's Binding.
 * @tparam Outer True to read the first warp slot, false for the second.
 * @tparam TrackPath Whether the stage accumulates path length, which the color
 *         stage consumes under HueMode::PATH_LENGTH.
 */
template <typename BindingT, bool Outer, bool TrackPath = false>
struct WarpProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static const auto &params(const FrameState &frame) {
    if constexpr (Outer)
      return frame.params.outer_warp;
    else
      return frame.params.inner_warp;
  }
  static auto prepare(const FrameState &frame) {
    using WarpT = std::remove_cvref_t<decltype(params(frame))>;
    if constexpr (std::is_same_v<WarpT, AffineParams>) {
      static_assert(
          requires { frame.params.source.lattice_cell_scale; },
          "the affine warp stage translates in lattice cells and requires a "
          "LatticeSourceParams source");
      return Pullback::Warp::prepare(
          params(frame), phase(frame), Outer ? frame.outer_rotation : 0.0f,
          1.0f / frame.params.source.lattice_cell_scale);
    } else if constexpr (std::is_same_v<WarpT, NoWarpParams> ||
                         std::is_same_v<WarpT, PolarParams>) {
      return Pullback::NoPrepared{};
    } else {
      return Pullback::Warp::prepare(params(frame), phase(frame));
    }
  }
  static float phase(const FrameState &frame) {
    if constexpr (Outer)
      return frame.outer_phase;
    else
      return frame.inner_phase;
  }
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.outer_noise;
  }
  static bool path_length_required(const FrameState &) { return TrackPath; }
};

/**
 * @brief Supplies the displacement field to the Pullback::Surface policies.
 * @tparam BindingT The effect's Binding.
 * @tparam TrackPath Whether the stage accumulates path length.
 * @pre The effect's `surface_type` is a surface-noise family, so the frame
 *      carries a noise pointer and a loop offset.
 */
template <typename BindingT, bool TrackPath = false> struct SurfaceProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.surface_noise;
  }
  static auto prepare(const FrameState &frame) {
    if constexpr (requires { frame.params.surface.direction; })
      return Pullback::Surface::prepare_direct(frame.surface_phase,
                                               frame.params.surface.direction);
    else
      return Pullback::Surface::prepare(frame.surface_phase);
  }
  static float scale(const FrameState &frame) {
    return frame.params.surface.scale;
  }
  static float strength(const FrameState &frame) {
    return frame.params.surface.strength;
  }
  static bool path_length_required(const FrameState &) { return TrackPath; }
};

/**
 * @brief Supplies the pattern and noise state to the Pullback::Source policies.
 * @details Covers every source family at once: the pattern accessors read the
 * prepared phases, the noise accessors the NoiseSourceParams fields. Only the
 * accessors an effect's chosen source policy names are instantiated.
 */
template <typename BindingT> struct SourceProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static const auto &params(const FrameState &frame) {
    return frame.params.source;
  }
  static Pullback::Source::PreparedSource prepare(const FrameState &frame) {
    return Pullback::Source::prepare(
        frame.source_primary, frame.source_secondary, frame.source_angle);
  }
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.source_noise;
  }
  static float noise_scale(const FrameState &frame) {
    return frame.params.source.noise_scale;
  }
  static float noise_time(const FrameState &frame) {
    return frame.source_noise_time;
  }
  static float noise_contrast(const FrameState &frame) {
    return frame.params.source.noise_contrast;
  }
};

/**
 * @brief Supplies the value-family fields to the Pullback::Transfer and
 *        Pullback::Coverage policies.
 * @details Names all three fields the value families define; an effect's material
 * stage instantiates only the accessors its transfer and coverage policies
 * call, so an IsoValueParams effect never touches `edge_width` and vice versa.
 */
template <typename BindingT> struct ValueProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static float iso_level(const FrameState &frame) {
    return frame.params.value.iso_level;
  }
  static float iso_width(const FrameState &frame) {
    return frame.params.value.iso_width;
  }
  static float edge_width(const FrameState &frame) {
    return frame.params.value.edge_width;
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
  static Pullback::Color::PaletteMappingWeights
  mapping_weights(const FrameState &frame) {
    return frame.palette_mapping;
  }
  static float mapping_frequency(const FrameState &frame) {
    return frame.params.color.mapping_frequency;
  }
  static float mapping_phase(const FrameState &frame) {
    return frame.params.color.mapping_phase;
  }
  static float oscillation_depth(const FrameState &frame) {
    return frame.params.color.phase_oscillation_depth;
  }
  static float oscillation_phase(const FrameState &frame) {
    return frame.palette_oscillation_phase;
  }
  static const BakedPalette &palette(const FrameState &frame) {
    return *frame.palette;
  }
  static Pullback::Color::HueMode hue_mode(const FrameState &) { return HueV; }
  static float hue_shift_amount(const FrameState &frame) {
    return frame.params.color.hue_shift_amount;
  }
  static Pullback::Color::HueRotationLutView
  hue_rotation(const FrameState &frame) {
    return {frame.hue_rotation_lut,
            hue_rotation_active<HueV>(frame.params.color)};
  }
  static Pullback::Color::HueNoiseLutView hue_noise(const FrameState &frame) {
    return {frame.hue_noise_lut,
            HueV == HueMode::NOISE &&
                hue_rotation_active<HueV>(frame.params.color)};
  }
  static Pullback::Color::BrightnessEnvelope
  brightness_envelope(const FrameState &) {
    return BrightnessV;
  }
  static float brightness_depth(const FrameState &frame) {
    return frame.params.color.brightness_depth;
  }
  static float opacity_low(const FrameState &frame) {
    return frame.params.color.opacity_low;
  }
  static float opacity_high(const FrameState &frame) {
    return frame.params.color.opacity_high;
  }
};

/**
 * @brief Interpolates one parameter family across a preset transition.
 * @details Driven by the family's field table; each field moves on the curve
 * its descriptor names. The Mobius coefficients and the palette mapping enum
 * snap to @p b once progress reaches 1.
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

inline ColorParams interpolate(const ColorParams &a, const ColorParams &b,
                               float t) {
  ColorParams value = Pullback::Fields::interpolate(a, b, t);
  value.palette_mapping = t < 1.0f ? a.palette_mapping : b.palette_mapping;
  return value;
}

/**
 * @brief Interpolates a whole parameter set, family by family.
 * @param from Parameters at progress 0.
 * @param to Parameters at progress 1.
 * @param progress Progress fraction, already eased by the caller.
 * @return The interpolated parameter set.
 */
template <typename SourceT, typename OuterWarpT, typename InnerWarpT,
          typename LensT, typename ValueT, typename SurfaceT>
inline Params<SourceT, OuterWarpT, InnerWarpT, LensT, ValueT, SurfaceT>
interpolate(
    const Params<SourceT, OuterWarpT, InnerWarpT, LensT, ValueT, SurfaceT>
        &from,
    const Params<SourceT, OuterWarpT, InnerWarpT, LensT, ValueT, SurfaceT> &to,
    float progress) {
  return {interpolate(from.source, to.source, progress),
          interpolate(from.projection, to.projection, progress),
          interpolate(from.outer_warp, to.outer_warp, progress),
          interpolate(from.inner_warp, to.inner_warp, progress),
          interpolate(from.surface, to.surface, progress),
          interpolate(from.lens, to.lens, progress),
          interpolate(from.value, to.value, progress),
          interpolate(from.color, to.color, progress)};
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
    if (!std::isfinite(value))
      return false;
  return true;
}

inline bool valid(const ColorParams &p) {
  return Pullback::Fields::valid(p) &&
         static_cast<uint8_t>(p.palette_mapping) <=
             static_cast<uint8_t>(Pullback::Color::PaletteMapping::REVERSE);
}

/**
 * @brief Whether every family of a parameter set is in range.
 * @return True only when all eight families pass.
 */
template <typename SourceT, typename OuterWarpT, typename InnerWarpT,
          typename LensT, typename ValueT, typename SurfaceT>
inline bool valid(const Params<SourceT, OuterWarpT, InnerWarpT, LensT, ValueT,
                               SurfaceT> &params) {
  return valid(params.source) && valid(params.projection) &&
         valid(params.outer_warp) && valid(params.inner_warp) &&
         valid(params.surface) && valid(params.lens) && valid(params.value) &&
         valid(params.color);
}

/** @brief Storage for one optional noise field; empty when disabled. */
template <bool Enabled> struct OptionalNoise {};
template <> struct OptionalNoise<true> {
  FastNoiseLite noise;
};

/** @brief Projection a composed effect's surface stage composes. */
enum class ProjectionKind : uint8_t {
  STEREOGRAPHIC,
  GNOMONIC_FOLDED,
  EQUIRECTANGULAR,
  FOLDED_SINUSOIDAL
};

/** @brief Transfer curve an effect's material stage composes. */
enum class TransferKind : uint8_t { LINEAR, ISO_CONTOUR };

/** @brief Coverage policy an effect's material stage composes. */
enum class CoverageKind : uint8_t { PROJECTION, PROJECTION_SQUARED, EDGE_FADE };

/**
 * @brief The pipeline choices an effect states beyond its parameter families.
 * @details Everything else about the six stages is derived: the source and
 * warp policies follow the parameter families, a Mobius lens family selects
 * the Mobius lens, a surface-noise family selects the curl displacement, and
 * path tracking follows HueMode::PATH_LENGTH.
 * @tparam ProjectionV Sphere-to-plane projection.
 * @tparam LensPolicyT Parameterless lens policy; ignored when the lens family
 *         carries Mobius coefficients.
 * @tparam TransferV Material transfer curve.
 * @tparam CoverageV Material coverage policy.
 */
template <ProjectionKind ProjectionV, typename LensPolicyT,
          TransferKind TransferV, CoverageKind CoverageV>
struct Spec {
  static constexpr ProjectionKind PROJECTION = ProjectionV;
  using LensPolicy = LensPolicyT;
  static constexpr TransferKind TRANSFER = TransferV;
  static constexpr CoverageKind COVERAGE = CoverageV;
};

template <typename Family, typename Binding> struct SourcePolicyFor;
template <typename B> struct SourcePolicyFor<GridSourceParams, B> {
  using Type = Pullback::Source::Grid<SourceProvider<B>>;
};
template <typename B> struct SourcePolicyFor<TwinWaveSourceParams, B> {
  using Type = Pullback::Source::TwinWave<SourceProvider<B>>;
};
template <typename B> struct SourcePolicyFor<SpiralSourceParams, B> {
  using Type = Pullback::Source::Spiral<SourceProvider<B>>;
};
template <typename B> struct SourcePolicyFor<LatticeSourceParams, B> {
  using Type = Pullback::Source::PrimitiveLattice<SourceProvider<B>>;
};

template <typename Family, typename Binding, bool Outer, bool TrackPath>
struct WarpPolicyFor;
template <typename B, bool O, bool T>
struct WarpPolicyFor<NoWarpParams, B, O, T> {
  using Type = Pullback::Warp::Identity;
};
template <typename B, bool O, bool T>
struct WarpPolicyFor<MirrorParams, B, O, T> {
  using Type = Pullback::Warp::MirrorTile<WarpProvider<B, O, T>>;
};
template <typename B, bool O, bool T>
struct WarpPolicyFor<WaveShearParams, B, O, T> {
  using Type = Pullback::Warp::WaveShear<WarpProvider<B, O, T>>;
};
template <typename B, bool O, bool T>
struct WarpPolicyFor<VectorNoiseParams, B, O, T> {
  using Type =
      Pullback::Warp::VectorNoise<WarpProvider<B, O, T>, NoiseBasis::SIMPLEX,
                                  Pullback::Warp::FlatEnvelope>;
};
template <typename B, bool O, bool T>
struct WarpPolicyFor<AffineParams, B, O, T> {
  using Type = Pullback::Warp::AffineFrame<WarpProvider<B, O, T>>;
};
template <typename B, bool O, bool T>
struct WarpPolicyFor<PolarParams, B, O, T> {
  using Type = Pullback::Warp::PolarChart<WarpProvider<B, O, T>,
                                          Pullback::Warp::LinearPolar, 1>;
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

template <typename LensFamily, typename SpecLens, typename Binding>
struct LensPolicyFor {
  using Type = SpecLens;
};
template <typename SpecLens, typename B>
struct LensPolicyFor<MobiusLensParams, SpecLens, B> {
  using Type = Pullback::Lens::Mobius<LensProvider<B>>;
};

template <TransferKind TransferV, typename Binding> struct TransferPolicyFor;
template <typename B> struct TransferPolicyFor<TransferKind::LINEAR, B> {
  using Type = Pullback::Transfer::Linear;
};
template <typename B> struct TransferPolicyFor<TransferKind::ISO_CONTOUR, B> {
  using Type = Pullback::Transfer::IsoContour<ValueProvider<B>>;
};

template <CoverageKind CoverageV, typename Binding> struct CoveragePolicyFor;
template <typename B> struct CoveragePolicyFor<CoverageKind::PROJECTION, B> {
  using Type = Pullback::Coverage::Projection;
};
template <typename B>
struct CoveragePolicyFor<CoverageKind::PROJECTION_SQUARED, B> {
  using Type = Pullback::Coverage::ProjectionSquared;
};
template <typename B> struct CoveragePolicyFor<CoverageKind::EDGE_FADE, B> {
  using Type = Pullback::Coverage::EdgeFade<ValueProvider<B>>;
};

/**
 * @brief A complete composed effect: the shared lifecycle plus a pipeline
 *        derived from the parameter families and a Spec.
 * @details The effect states its families, its Spec and its identity
 * constants; every stage typedef, the render pipeline, shade() and the shared
 * lifecycle — parameter registration, preset choreography, palette cycling,
 * camera walks and noise clocks — are assembled here. Required `Derived`
 * members are `PRESET_IDS`, `PARAMETER_SCHEMA_VERSION`, `PRESET_DWELL_FRAMES`,
 * `ANIMATED_PROJECTION`, `initial_params` and `preset_params`, plus an
 * `OUTER_NOISE_SEED` / `SOURCE_NOISE_SEED` / `SURFACE_NOISE_SEED` for each
 * noise field the parameter set and `Has*Noise` flags request. Optional
 * members, detected by `requires` and defaulted when absent, are
 * `ANIMATED_MOBIUS`, `USES_CENTRAL_MERIDIAN` and an `after_composed_init()`
 * hook. An effect wanting a different shade emission shadows shade() with the
 * same body under its own attribute. A surface-noise family places its
 * surface stage out of line in flash.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam Derived The effect class deriving from this base.
 * @tparam ParamsT The effect's `Pullback::Params` specialization.
 * @tparam SpecT The effect's `Spec`.
 * @tparam Harmony Palette harmony the generated palettes are drawn from.
 * @tparam HueV Hue-rotation source: none, noise field, or path length.
 * @tparam BrightnessV Brightness envelope applied by the color stage.
 * @tparam HasOuterNoise Whether the outer-camera stage owns a noise field.
 * @tparam HasSourceNoise Whether the source stage owns a noise field.
 */
template <int W, int H, typename Derived, typename ParamsT, typename SpecT,
          PaletteHarmony Harmony, HueMode HueV,
          Pullback::Color::BrightnessEnvelope BrightnessV,
          bool HasOuterNoise = false, bool HasSourceNoise = false>
class ComposedEffect : public ChoreographedEffect<Derived, ParamsT> {
  using Choreography = ChoreographedEffect<Derived, ParamsT>;
  friend Choreography;

public:
  using Params = ParamsT;
  using Frame = Pullback::FrameState<ParamsT>;
  using Binding = Pullback::Binding<Frame>;
  /** Preset policy: an automatic change crossfades the parameters; pause
      never freezes an in-flight crossfade. */
  static constexpr Segue::Lerp PRESET_SEGUE{480, ease_in_out_sin};
  /** Whether the effect displaces its surface, and so owns a surface noise
      field and needs a `Derived::SURFACE_NOISE_SEED`. */
  static constexpr bool HAS_SURFACE_NOISE =
      !std::is_same_v<typename ParamsT::surface_type, NoSurfaceParams>;

private:
  static constexpr bool TRACK_PATH = HueV == HueMode::PATH_LENGTH;
  static constexpr bool DIRECT_SURFACE =
      std::is_same_v<typename ParamsT::surface_type, DirectSurfaceParams>;
  using SurfacePolicy =
      std::conditional_t<HAS_SURFACE_NOISE && !DIRECT_SURFACE,
                         Pullback::Surface::CurlNoise<
                             SurfaceProvider<Binding, TRACK_PATH>,
                             NoiseBasis::SIMPLEX, Pullback::Surface::Euler>,
                         Pullback::Surface::Identity>;
  using PostLensSurfacePolicy = std::conditional_t<
      DIRECT_SURFACE,
      Pullback::Surface::DirectNoise<SurfaceProvider<Binding, TRACK_PATH>,
                                     NoiseBasis::SIMPLEX>,
      Pullback::Surface::Identity>;
  using SurfaceImplementation = Pullback::Stage::SurfaceProject<
      Binding, SurfacePolicy,
      typename LensPolicyFor<typename ParamsT::lens_type,
                             typename SpecT::LensPolicy, Binding>::Type,
      PostLensSurfacePolicy,
      typename ProjectionPolicyFor<SpecT::PROJECTION, Binding>::Type>;

public:
  /// Pullback stages, ordered from the view vector back to the source field.
  using OuterCameraStage =
      Pullback::Stage::OuterCamera<Binding, OuterCameraProvider<Binding>>;
  using SurfaceStage = std::conditional_t<
      HAS_SURFACE_NOISE,
      Pullback::Stage::Placed<Pullback::CodeEmission::OUT_OF_LINE_FLASH,
                              SurfaceImplementation>,
      SurfaceImplementation>;
  using PlanarWarpStage = Pullback::Stage::PlanarWarp<
      Binding,
      typename WarpPolicyFor<typename ParamsT::outer_warp_type, Binding, true,
                             TRACK_PATH>::Type,
      typename WarpPolicyFor<typename ParamsT::inner_warp_type, Binding, false,
                             TRACK_PATH>::Type>;
  using SourceStage = Pullback::Stage::Source<
      Binding,
      typename SourcePolicyFor<typename ParamsT::source_type, Binding>::Type>;
  using MaterialStage = Pullback::Stage::Material<
      Binding, Pullback::Weight::Projection,
      typename TransferPolicyFor<SpecT::TRANSFER, Binding>::Type,
      typename CoveragePolicyFor<SpecT::COVERAGE, Binding>::Type>;
  using ColorStage =
      Pullback::Stage::Color<Binding,
                             Pullback::Color::GeneratedPalette<
                                 ColorProvider<Binding, HueV, BrightnessV>>>;
  using RenderPipeline =
      Pullback::Pipeline<Binding, OuterCameraStage, SurfaceStage,
                         PlanarWarpStage, SourceStage, MaterialStage,
                         ColorStage>;
  using FrameState = typename RenderPipeline::Frame;

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
    configure_presets(Derived::PRESET_IDS.size());
    state = persistent_arena.make<State>();
    use_parameter_storage(persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);
    configure_noise(state->color_noise, 6047);
    if constexpr (HasOuterNoise)
      configure_noise(state->outer.noise, Derived::OUTER_NOISE_SEED);
    if constexpr (HasSourceNoise)
      configure_noise(state->source.noise, Derived::SOURCE_NOISE_SEED);
    if constexpr (HAS_SURFACE_NOISE)
      configure_noise(state->surface.noise, Derived::SURFACE_NOISE_SEED);
    init_gamut_lut(persistent_arena, GAMUT_LUT_ANGLE_STEPS, GAMUT_LUT_L_STEPS);
    palette_cycler.init_generated(persistent_arena, next_palette, this, 0, 600,
                                  ease_in_out_sin);
    timeline.add(0, Animation::RandomWalk<W>(projection_walk, UP,
                                             state->projection_walk_noise));
    timeline.add(
        0, Animation::RandomWalk<W>(outer_walk, UP, state->outer_walk_noise));
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
      begin_automatic_transition();
      advance_runtime();
      update_spatial_frames();
      update_palette_chroma();
      palette_cycler.step();
    }
    const typename Derived::RenderPipeline::Frame frame =
        Derived::RenderPipeline::prepare(prepare_frame());
    {
      HS_PROFILE(fx_shader_draw);
      Scan::Shader::draw_cached<W, H, 1>(canvas, [&frame](const Vector &view) {
        return Derived::shade(view, frame);
      });
    }
  }

  /**
   * @brief Shades one pixel through the fully inlined pipeline.
   * @param view Unit view direction for the pixel.
   * @param frame Per-frame transforms, params and LUTs from the runtime.
   */
  static HS_FLASH_INLINE Color4 shade(const Vector &view,
                                      const FrameState &frame) {
    return RenderPipeline::shade(view, frame);
  }

  /** @brief Whether a parameter set is admissible, family by family. */
  static bool valid_params(const Params &params) {
    return Pullback::valid(params);
  }

protected:
  /** Descriptors the arena-backed parameter array holds; every slider an effect
      registers has to fit. */
  static constexpr size_t PARAM_CAPACITY = 48;

  using Choreography::anims_paused;
  using Choreography::begin_automatic_transition;
  using Choreography::configure_presets;
  using Choreography::params;
  using Choreography::register_animated_param;
  using Choreography::register_fields;
  using Choreography::timeline;
  using Choreography::transition;
  using Choreography::use_parameter_storage;

  /** @brief Adopts a snap target and re-derives the palette mapping weights. */
  HS_COLD_MEMBER void adopt_params(const Params &target) {
    params = target;
    palette_mapping = Pullback::Color::PaletteMappingWeights::single(
        target.color.palette_mapping);
  }

  /** @brief Captures the palette-mapping endpoints of an arming crossfade. */
  HS_COLD_MEMBER void transition_armed(const Params &target) {
    mapping_from = palette_mapping;
    mapping_to = Pullback::Color::PaletteMappingWeights::single(
        target.color.palette_mapping);
  }

  /**
   * @brief Writes the interpolated parameters of an in-flight transition.
   * @details Under `Derived::ANIMATED_MOBIUS` the live lens coefficients are
   * carried across the interpolation, so the transition does not overwrite the
   * timeline animation driving them.
   */
  HS_COLD_MEMBER void blend_params(float progress) {
    MobiusParams animated_mobius;
    if constexpr (requires { Derived::ANIMATED_MOBIUS; })
      if constexpr (Derived::ANIMATED_MOBIUS)
        animated_mobius = params.lens.mobius;
    params = Pullback::interpolate(transition.from, transition.to, progress);
    if constexpr (requires { Derived::ANIMATED_MOBIUS; })
      if constexpr (Derived::ANIMATED_MOBIUS)
        params.lens.mobius = animated_mobius;
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
  HS_COLD_MEMBER void start_mobius_animation(float scale, int duration) {
    if constexpr (requires { params.lens.mobius; })
      timeline.add_pausable(0,
                            Animation::MobiusWarpCircular(
                                params.lens.mobius, scale, duration, true),
                            &anims_paused);
  }

private:
  struct State {
    std::array<Pixel, Pullback::Color::HueRotationLutView::SIZE>
        hue_rotation_lut;
    std::array<int8_t, Pullback::Color::HueNoiseLutView::SIZE> hue_noise_lut;
    FastNoiseLite color_noise;
    OptionalNoise<HasOuterNoise> outer;
    OptionalNoise<HasSourceNoise> source;
    OptionalNoise<HAS_SURFACE_NOISE> surface;
    FastNoiseLite projection_walk_noise;
    FastNoiseLite outer_walk_noise;
    /** Inputs the hue-noise LUT was last baked from; the negative seeds match
        no admissible parameter, forcing the first bake. */
    float hue_noise_lut_scale = -1.0f;
    float hue_noise_lut_phase = -1.0f;
  };

  // init() takes the gamut LUT, the palette cycler's generated arena, the
  // external ParamDef array and one State from the persistent arena.
  static constexpr size_t FOOTPRINT_BYTES =
      gamut_lut_bytes(GAMUT_LUT_ANGLE_STEPS, GAMUT_LUT_L_STEPS) +
      PaletteCycler::generated_arena_bytes() +
      PARAM_CAPACITY * sizeof(ParamDef) + sizeof(State) + alignof(State);
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
      "Pullback::ComposedEffect persistent footprint exceeds the default "
      "partition");

  static void configure_noise(FastNoiseLite &noise, int32_t seed) {
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise.SetSeed(seed);
    noise.SetFrequency(1.0f);
  }

  /** @brief Whether a gated field's slider exists for this effect. */
  static constexpr bool field_gate_open(Pullback::FieldGate gate) {
    switch (gate) {
    case Pullback::FieldGate::ALWAYS:
      return true;
    case Pullback::FieldGate::ANIMATED_PROJECTION:
      return Derived::ANIMATED_PROJECTION;
    case Pullback::FieldGate::CENTRAL_MERIDIAN:
      if constexpr (requires { Derived::USES_CENTRAL_MERIDIAN; })
        return Derived::USES_CENTRAL_MERIDIAN;
      else
        return false;
    }
    return true;
  }

  /**
   * @brief Registers one warp slot: the slot-named speed slider, then the
   *        family's named fields.
   */
  template <typename T>
  HS_COLD_MEMBER void register_warp_fields(T &warp, const char *slot_name) {
    if constexpr (!std::is_same_v<T, NoWarpParams>) {
      static_assert(T::FIELDS[0].member == &T::speed,
                    "warp slot registration expects speed first");
      register_animated_param(slot_name, &warp.speed, T::FIELDS[0].min,
                              T::FIELDS[0].max);
      register_fields(warp);
    }
  }

  HS_COLD_MEMBER void register_parameters() {
    register_fields(params.source);
    register_fields(params.projection);
    register_fields(params.surface);
    register_warp_fields(params.outer_warp, "Planar Warp 1 Speed");
    register_warp_fields(params.inner_warp, "Planar Warp 2 Speed");
    register_fields(params.value);
    if constexpr (requires { params.lens.mobius; }) {
      register_animated_param("Mobius A Re", &params.lens.mobius.a.re, -4.0f,
                              4.0f);
      register_animated_param("Mobius A Im", &params.lens.mobius.a.im, -4.0f,
                              4.0f);
      register_animated_param("Mobius B Re", &params.lens.mobius.b.re, -4.0f,
                              4.0f);
      register_animated_param("Mobius B Im", &params.lens.mobius.b.im, -4.0f,
                              4.0f);
      register_animated_param("Mobius C Re", &params.lens.mobius.c.re, -4.0f,
                              4.0f);
      register_animated_param("Mobius C Im", &params.lens.mobius.c.im, -4.0f,
                              4.0f);
      register_animated_param("Mobius D Re", &params.lens.mobius.d.re, -4.0f,
                              4.0f);
      register_animated_param("Mobius D Im", &params.lens.mobius.d.im, -4.0f,
                              4.0f);
    }
    register_animated_param("Palette Chroma", &params.color.palette_chroma,
                            0.0f, 1.0f);
    register_animated_param("Palette Mapping", &params.color.palette_mapping,
                            PALETTE_MAPPING_OPTIONS,
                            PALETTE_MAPPING_EXPORT_OPTIONS,
                            std::size(PALETTE_MAPPING_OPTIONS));
    register_animated_param("Mapping Frequency",
                            &params.color.mapping_frequency, 1.0f, 32.0f);
    register_animated_param("Mapping Phase", &params.color.mapping_phase, -1.0f,
                            1.0f);
    register_animated_param("Phase Oscillation Depth",
                            &params.color.phase_oscillation_depth, 0.0f, 1.0f);
    register_animated_param("Phase Oscillation Speed",
                            &params.color.phase_oscillation_speed, -0.01f,
                            0.01f);
    if constexpr (BrightnessV != Pullback::Color::BrightnessEnvelope::NONE)
      register_animated_param("Brightness Depth",
                              &params.color.brightness_depth, 0.0f, 1.0f);
    register_animated_param("Value Opacity Low", &params.color.opacity_low,
                            0.0f, 1.0f);
    register_animated_param("Value Opacity High", &params.color.opacity_high,
                            0.0f, 1.0f);
    register_animated_param("Hue Shift Amount", &params.color.hue_shift_amount,
                            HueV == HueMode::PATH_LENGTH ? -4.0f : -1.0f,
                            HueV == HueMode::PATH_LENGTH ? 4.0f : 1.0f);
    if constexpr (HueV == HueMode::NOISE) {
      register_animated_param("Hue Noise Scale", &params.color.hue_noise_scale,
                              1.0f / 64.0f, 8.0f);
      register_animated_param("Hue Noise Speed", &params.color.hue_noise_speed,
                              -0.001f, 0.001f);
    }
  }

  /**
   * @brief Steps every phase clock the effect's parameter families define.
   * @details Each clock is compiled in only when its field exists, so an effect
   * pays for exactly the phases its stages read. The affine frame rotation is
   * accumulated for the outer warp slot alone.
   */
  HS_COLD_MEMBER void advance_runtime() {
    if constexpr (requires { params.source.speed; }) {
      source_primary = fmodf(source_primary + params.source.speed, TWO_PI_F);
      if constexpr (requires { params.source.secondary_rate; })
        source_secondary =
            fmodf(source_secondary +
                      params.source.speed * params.source.secondary_rate,
                  TWO_PI_F);
      if constexpr (requires { params.source.angle_rate; })
        source_angle = fmodf(source_angle + params.source.angle_rate, TWO_PI_F);
    }
    if constexpr (requires { params.source.noise_time_rate; })
      source_noise_time =
          wrap_t(source_noise_time + params.source.noise_time_rate);
    if constexpr (HAS_SURFACE_NOISE)
      surface_phase = wrap_t(surface_phase + params.surface.speed);
    if constexpr (Derived::ANIMATED_PROJECTION)
      projection_spin =
          fmodf(projection_spin + params.projection.spin_rate, TWO_PI_F);
    hue_noise_phase = wrap_t(hue_noise_phase + params.color.hue_noise_speed);
    if constexpr (std::is_same_v<typename Params::outer_warp_type,
                                 AffineParams>)
      outer_rotation =
          TWO_PI_F *
          wrap_t((outer_rotation +
                  params.outer_warp.speed * params.outer_warp.rotation_rate) /
                 TWO_PI_F);
    outer_phase = wrap_t(outer_phase + params.outer_warp.speed);
    inner_phase = wrap_t(inner_phase + params.inner_warp.speed);
    palette_oscillation_phase = wrap_t(palette_oscillation_phase +
                                       params.color.phase_oscillation_speed);
  }

  HS_COLD_MEMBER void update_spatial_frames() {
    // prepare_frame() reads projection_conjugate only for an animated
    // projection.
    if constexpr (Derived::ANIMATED_PROJECTION) {
      const Quaternion projection = projection_walk.get();
      const Quaternion projection_delta =
          projection * projection_walk_previous.conjugate();
      projection_walk_previous = projection;
      projection_wander = (scaled_rotation_delta(projection_delta.normalized(),
                                                 params.projection.wander) *
                           projection_wander)
                              .normalized();
      projection_conjugate = (make_rotation(Y_AXIS, projection_spin) *
                              base_orientation * projection_wander)
                                 .conjugate();
    }
    const Quaternion outer = outer_walk.get();
    const Quaternion outer_delta = outer * outer_walk_previous.conjugate();
    outer_walk_previous = outer;
    outer_wander = (scaled_rotation_delta(outer_delta.normalized(),
                                          params.projection.camera_wander) *
                    outer_wander)
                       .normalized();
    outer_conjugate = outer_wander.conjugate();
  }

  /**
   * @brief Bakes the frame's LUTs and snapshots everything the scan reads.
   * @details The hue-noise LUT is rebuilt only when its scale or phase moved,
   * while the hue-rotation LUT is rebuilt on exactly the condition
   * hue_rotation_active() reports, which is also the flag the returned frame
   * hands the color stage.
   * @return The frame state for this draw, valid until the next draw_frame().
   */
  HS_COLD_MEMBER Frame prepare_frame() {
    HS_PROFILE(fx_prepare_frame);
    if constexpr (HueV == HueMode::NOISE) {
      if (state->hue_noise_lut_scale != params.color.hue_noise_scale ||
          state->hue_noise_lut_phase != hue_noise_phase) {
        Pullback::Color::prepare_hue_noise_lut(
            std::span<int8_t, Pullback::Color::HueNoiseLutView::SIZE>(
                state->hue_noise_lut),
            state->color_noise, params.color.hue_noise_scale, hue_noise_phase);
        state->hue_noise_lut_scale = params.color.hue_noise_scale;
        state->hue_noise_lut_phase = hue_noise_phase;
      }
    }
    if (hue_rotation_active<HueV>(params.color))
      Pullback::Color::prepare_hue_rotation_lut(
          std::span<Pixel, Pullback::Color::HueRotationLutView::SIZE>(
              state->hue_rotation_lut),
          palette_cycler.palette());
    const FastNoiseLite *outer_noise = nullptr;
    const FastNoiseLite *source_noise = nullptr;
    if constexpr (HasOuterNoise)
      outer_noise = &state->outer.noise;
    if constexpr (HasSourceNoise)
      source_noise = &state->source.noise;
    const FastNoiseLite *surface_noise = nullptr;
    if constexpr (HAS_SURFACE_NOISE)
      surface_noise = &state->surface.noise;
    return {Derived::ANIMATED_PROJECTION ? projection_conjugate : Quaternion(),
            outer_conjugate,
            outer_noise,
            source_noise,
            surface_noise,
            &palette_cycler.palette(),
            state->hue_rotation_lut.data(),
            state->hue_noise_lut.data(),
            params,
            palette_mapping,
            source_primary,
            source_secondary,
            source_angle,
            outer_phase,
            inner_phase,
            outer_rotation,
            source_noise_time,
            surface_phase,
            palette_oscillation_phase};
  }

  HS_COLD_MEMBER void update_palette_chroma() {
    if (palette_chroma == params.color.palette_chroma)
      return;
    palette_chroma = params.color.palette_chroma;
    palette_cycler.set_generated_chroma(palette_chroma);
  }

  /**
   * @brief Generator the palette cycler calls for each palette in the cycle.
   * @details Advances the hue on every palette after the first, so the opening
   * palette is the one the effect's authored hue names.
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
        effect.params.color.palette_chroma)};
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
          params.color.palette_mapping);
  Pullback::Color::PaletteMappingWeights mapping_from;
  Pullback::Color::PaletteMappingWeights mapping_to;
  Orientation<> projection_walk;
  Orientation<> outer_walk;
  Quaternion projection_walk_previous;
  Quaternion outer_walk_previous;
  Quaternion projection_wander;
  Quaternion outer_wander;
  Quaternion projection_conjugate;
  Quaternion outer_conjugate;
  Quaternion base_orientation =
      make_rotation(Vector(0, 0, -1), Vector(0, -1, 0));
  float source_primary = 0.0f;
  float source_secondary = 0.0f;
  float source_angle = 0.0f;
  float source_noise_time = 0.0f;
  float surface_phase = 0.0f;
  float projection_spin = 0.0f;
  float hue_noise_phase = 0.0f;
  float outer_phase = 0.0f;
  float inner_phase = 0.0f;
  float outer_rotation = 0.0f;
  float palette_oscillation_phase = 0.0f;
  float palette_chroma = -1.0f;
  uint32_t palette_hue = 0;
  PaletteCycler palette_cycler;
};

} // namespace Pullback
