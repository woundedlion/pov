/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file FixedLookRuntime.h
 * @brief Shared machinery for the fixed-look effect family: the parameter
 *        families a look composes into a `Params`, the frame-state providers
 *        its pullback pipeline reads, and the `Runtime` base that drives them.
 */

#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <type_traits>

#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"
#include "core/engine/fixed_pipeline.h"
#include "core/math/noise_field.h"
#include "core/render/pullback.h"

namespace FixedLook {

/** @brief What drives the color stage's hue rotation, if anything. */
enum class HueMode : uint8_t {
  NONE,       /**< No hue rotation; the palette color is used as sampled. */
  NOISE,      /**< Rotation amount read from a cube-face noise LUT. */
  PATH_LENGTH /**< Rotation amount read from the accumulated path length. */
};

/**
 * @brief Source parameters for the coupled sine grid
 *        (Pullback::Source::Grid).
 */
struct GridSourceParams {
  float pattern_freq = 1.0f; /**< Scale applied to the warped plane
                                   coordinates before the grid is sampled. */
  float speed = 0.0f;        /**< Per-frame advance of the primary phase. */
  float complexity = 0.0f;   /**< Amount of cross-axis coupling folded into the
                                   grid coordinates. */
  float pattern_mix = 0.0f;  /**< Blend from the coupled pattern at 0 to the
                                   direct sine product at 1. */
  float secondary_rate = 0.0f; /**< Secondary phase rate, as a multiple of
                                    `speed`. */
  float angle_rate = 0.0f;     /**< Per-frame advance of the source rotation. */
};

/**
 * @brief Source parameters for the two-wave interference field
 *        (Pullback::Source::TwinWave).
 */
struct TwinWaveSourceParams {
  float pattern_freq = 1.0f;   /**< Plane-coordinate scale before sampling. */
  float speed = 0.0f;          /**< Per-frame advance of the primary phase. */
  float secondary_rate = 0.0f; /**< Secondary phase rate, as a multiple of
                                    `speed`. */
  float angle_rate = 0.0f; /**< Per-frame advance of the angle between the two
                                waves. */
};

/**
 * @brief Source parameters for the noise-contour sources
 *        (Pullback::Source::ProjectedNoise and SphericalNoise).
 */
struct NoiseSourceParams {
  float noise_scale = 1.0f;    /**< Spatial scale of the sampled field. */
  float noise_contrast = 0.0f; /**< Contour sharpening applied to the sample. */
  float noise_time_rate = 0.0f; /**< Per-frame advance of the noise time
                                     coordinate. */
};

/**
 * @brief Source parameters for the per-cell primitive lattice
 *        (Pullback::Source::PrimitiveLattice).
 */
struct LatticeSourceParams {
  float lattice_cell_scale = 1.0f;  /**< Lattice cells per plane unit. */
  float lattice_shape_blend = 0.0f; /**< Cell primitive, from a circle at 0 to a
                                         rounded square at 1. */
  float lattice_softness = 0.05f;   /**< Half-width of the ramp across the
                                         primitive's boundary. */
  float lattice_radius = 0.25f;     /**< Primitive radius in cell units. */
};

/** @brief Surface stage placeholder for a look that carries no displacement. */
struct NoSurfaceParams {};

/**
 * @brief Surface parameters for the sphere-space noise displacements
 *        (Pullback::Surface::CurlNoise and DirectNoise).
 */
struct SurfaceNoiseParams {
  float scale = 1.0f;    /**< Spatial scale of the displacement field. */
  float strength = 0.0f; /**< Displacement distance; 0 skips the stage. */
  float speed = 0.0f;    /**< Per-frame advance of the field's loop phase. */
};

/** @brief Projection and camera parameters, shared by every look. */
struct ProjectionParams {
  float pole_fade = 1.0f; /**< Falloff applied to the projection's radial
                                attenuation. */
  float spin_rate = 0.0f; /**< Per-frame spin of the projection frame about Y;
                                only read under `ANIMATED_PROJECTION`. */
  float wander = 0.0f;    /**< Fraction of the projection random-walk delta
                                absorbed each frame. */
  float camera_wander = 0.0f;    /**< Same, for the outer camera random walk. */
  float central_meridian = 0.0f; /**< Central meridian handed to projections
                                      that take one, in radians. */
};

/**
 * @brief Warp slot placeholder for a look whose warp policy is an identity.
 * @details `speed` still advances that slot's phase clock, so a look can drive
 * a phase it exposes no warp for.
 */
struct NoWarpParams {
  float speed = 0.0f; /**< Per-frame advance of the slot's phase. */
};

/**
 * @brief Warp parameters for the mirrored tiling
 *        (Pullback::Warp::MirrorTile).
 */
struct MirrorParams {
  float speed = 0.0f;    /**< Per-frame advance of the slot's phase. */
  float rotation = 0.0f; /**< Rotation of the fold lattice, in radians. */
  float cell_x = 1.0f;   /**< Mirror cell width in plane units. */
  float cell_y = 1.0f;   /**< Mirror cell height in plane units. */
  float offset_x = 0.0f; /**< Pre-fold translation along x; scrolls with the
                              slot's phase. */
  float offset_y = 0.0f; /**< Pre-fold translation along y; does not scroll. */
};

/** @brief Warp parameters for the sine shear (Pullback::Warp::WaveShear). */
struct WaveShearParams {
  float speed = 0.0f;       /**< Per-frame advance of the slot's phase. */
  float strength = 0.0f;    /**< Shear amplitude; 0 skips the stage. */
  float frequency = 1.0f;   /**< Spatial frequency along the field axis. */
  float field_angle = 0.0f; /**< Field axis direction, in radians. */
  float edge_width = 0.1f;  /**< Fade band width, read only under an
                                 EdgeFadeEnvelope. */
};

/**
 * @brief Warp parameters for the noise-vector displacement
 *        (Pullback::Warp::VectorNoise).
 */
struct VectorNoiseParams {
  float speed = 0.0f;        /**< Per-frame advance of the slot's phase, which
                                  walks the noise loop. */
  float strength = 0.0f;     /**< Displacement amplitude; 0 skips the stage. */
  float scale = 1.0f;        /**< Spatial scale of the sampled field. */
  float vector_angle = 0.0f; /**< Rotation applied to the sampled vector, in
                                  radians. */
  float edge_width = 0.1f;   /**< Fade band width, read only under an
                                  EdgeFadeEnvelope. */
};

/**
 * @brief Warp parameters for the affine frame change
 *        (Pullback::Warp::AffineFrame).
 * @details Translation is expressed in lattice cells, so this family is only
 * valid alongside a LatticeSourceParams source.
 */
struct AffineParams {
  float speed = 0.0f;         /**< Per-frame advance of the slot's phase. */
  float rotation_rate = 0.0f; /**< Frame rotation rate; read only in the outer
                                   slot. */
  float translation_x = 0.0f; /**< Translation along x, in lattice cells per
                                   phase turn. */
  float translation_y = 0.0f; /**< Translation along y, in lattice cells per
                                   phase turn. */
  float scale_x = 1.0f; /**< Scale along x, oscillated over the phase cycle. */
  float scale_y = 1.0f; /**< Scale along y, oscillated over the phase cycle. */
  float shear = 0.0f;   /**< Shear, oscillated over the phase cycle. */
};

/** @brief Warp parameters for the polar chart (Pullback::Warp::PolarChart). */
struct PolarParams {
  float speed = 0.0f;         /**< Per-frame advance of the slot's phase, which
                                   offsets the angular coordinate. */
  float radial_scale = 1.0f;  /**< Scale applied to the radial coordinate. */
  float radial_phase = 0.0f;  /**< Offset added to the radial coordinate. */
  float angular_phase = 0.0f; /**< Offset added to the angular coordinate. */
};

/** @brief Lens placeholder for a look with no parameterized lens. */
struct NoLensParams {};
/** @brief Lens parameters for the Mobius map (Pullback::Lens::Mobius). */
struct MobiusLensParams {
  /** Mobius coefficients; the default is the identity map. */
  MobiusParams mobius{0.7071067811865475f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                      0.7071067811865475f, 0.0f};
};

/** @brief Value placeholder for a material stage that takes no parameters. */
struct LinearValueParams {};
/** @brief Value parameters for the edge fade (Pullback::Coverage::EdgeFade). */
struct EdgeValueParams {
  /** Fade band width in the projection's edge-distance units; 0 makes the edge
      a hard cut. */
  float edge_width = 0.1f;
};
/**
 * @brief Value parameters for the iso band
 *        (Pullback::Transfer::IsoContour).
 */
struct IsoValueParams {
  float iso_level = 0.5f;  /**< Source value the band is centered on. */
  float iso_width = 0.05f; /**< Half-width of the band's plateau. */
};

/** @brief Palette and hue parameters, shared by every look. */
struct ColorParams {
  float hue_shift_amount = 0.0f; /**< Hue rotation magnitude; 0 disables the
                                      rotation entirely. */
  float hue_noise_scale = 1.0f;  /**< Spatial scale of the hue-noise LUT. */
  float hue_noise_speed = 0.0f;  /**< Per-frame advance of the hue-noise loop
                                      phase; a change rebuilds the LUT. */
  float palette_chroma = 0.62f;  /**< Chroma the generated palettes are baked
                                      at. */
  /** Palette repeats across the value range. */
  float mapping_frequency = 1.0f;
  float mapping_phase = 0.0f;           /**< Offset into the palette. */
  float phase_oscillation_depth = 0.0f; /**< Amplitude of the sinusoidal wobble
                                             added to `mapping_phase`. */
  /** Per-frame advance of that wobble. */
  float phase_oscillation_speed = 0.0f;
  float brightness_depth = 1.0f; /**< Depth of the brightness envelope; only
                                      registered when the envelope is not
                                      NONE. */
  float opacity_low = 1.0f;      /**< Alpha gain at source value 0. */
  float opacity_high = 1.0f;     /**< Alpha gain at source value 1. */
  /** Palette mapping curve; snapped, not blended, by interpolate(). */
  Pullback::Color::PaletteMapping palette_mapping =
      Pullback::Color::PaletteMapping::LINEAR;
};

/**
 * @brief One look's complete parameter set, one family per pipeline stage.
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

/** @brief The source stage's phases, resolved once per frame. */
struct PreparedSource {
  float primary;   /**< Primary phase, wrapped into [0,2pi). */
  float secondary; /**< Secondary phase, wrapped into [0,2pi). */
  float angle;     /**< Source rotation, in radians. */
  float angle_cos; /**< Cosine of `angle`. */
  float angle_sin; /**< Sine of `angle`. */
};

/** @brief Affine warp coefficients, with the phase oscillation applied. */
struct PreparedAffine {
  float translation_x; /**< Translation along x, in plane units. */
  float translation_y; /**< Translation along y, in plane units. */
  float scale_x;       /**< Scale along x at this frame's phase. */
  float scale_y;       /**< Scale along y at this frame's phase. */
  float shear;         /**< Shear at this frame's phase. */
};

/** @brief Mirror warp offsets, with the phase scroll already folded in. */
struct PreparedMirror {
  float offset_x; /**< Pre-fold translation along x. */
  float offset_y; /**< Pre-fold translation along y. */
};

/** @brief This frame's point on the noise field's closed loop. */
struct PreparedNoiseLoop {
  float diagonal; /**< Offset added to both planar noise coordinates. */
  float z;        /**< Third noise coordinate. */
};

/**
 * @brief The transform half of PreparedWarp, one alternative per warp family.
 * @details The active member follows the slot's parameter family: MirrorParams
 * uses `mirror`, VectorNoiseParams `noise_loop`, AffineParams `affine`. The
 * remaining families leave the union zeroed, and their warp policies do not
 * read it.
 */
union PreparedWarpTransform {
  PreparedAffine affine;
  PreparedMirror mirror;
  PreparedNoiseLoop noise_loop;
};

/** @brief One planar warp slot's per-frame state. */
struct PreparedWarp {
  float rotation_cos;              /**< Cosine of the slot's rotation angle. */
  float rotation_sin;              /**< Sine of the slot's rotation angle. */
  PreparedWarpTransform transform; /**< Family-specific coefficients. */
};

/** @brief Surface stage per-frame state; empty unless the look displaces. */
template <typename SurfaceT> struct PreparedSurface {};
template <> struct PreparedSurface<SurfaceNoiseParams> {
  const FastNoiseLite *noise; /**< The runtime's surface noise field. */
  Vector loop_offset;         /**< This frame's point on the field's loop. */
};

/**
 * @brief Everything one frame's shading reads, resolved before the scan.
 * @details Handed to `Derived::shade` by value and read through the providers
 * below. Its pointers alias the runtime's persistent state and the palette
 * cycler's current bake, so a frame outlives only the draw_frame() call that
 * built it.
 * @tparam ParamsT The look's `FixedLook::Params` specialization.
 */
template <typename ParamsT> struct FrameState {
  /** Conjugate of the projection orientation; identity unless the look sets
      `ANIMATED_PROJECTION`. */
  Quaternion projection_conjugate;
  /** Conjugate of the outer camera orientation. */
  Quaternion outer_conjugate;
  PreparedSource prepared_source;
  PreparedWarp prepared_outer;       /**< State for the first warp slot. */
  PreparedWarp prepared_inner;       /**< State for the second warp slot. */
  const FastNoiseLite *outer_noise;  /**< Null unless `HasOuterNoise`. */
  const FastNoiseLite *source_noise; /**< Null unless `HasSourceNoise`. */
  const BakedPalette *palette;       /**< The cycler's current bake. */
  /** Hue-rotation LUT base; current only when hue_rotation_active(). */
  const Pixel *hue_rotation_lut;
  /** Hue-noise LUT base; current only under HueMode::NOISE with an active
      rotation. */
  const int8_t *hue_noise_lut;
  ParamsT params; /**< The frame's parameter values, already interpolated if a
                       preset transition is in flight. */
  /** Palette mapping weights, blended across a preset transition. */
  Pullback::Color::PaletteMappingWeights palette_mapping;
  float outer_phase;               /**< First warp slot's phase. */
  float inner_phase;               /**< Second warp slot's phase. */
  float source_noise_time;         /**< Source noise time coordinate. */
  float palette_oscillation_phase; /**< Phase of the mapping wobble. */
  PreparedSurface<typename ParamsT::surface_type> surface;
};

/**
 * @brief Ties a pullback pipeline to one look's frame state.
 * @details Fixed looks render uninstrumented, so the binding pins
 * Pullback::NoInstrumentation.
 * @tparam FrameT The look's FrameState specialization.
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
 * @details Exposes every accessor the projection policies name; a look pays
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
 * @tparam BindingT The look's Binding.
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
  static const PreparedWarp &prepared(const FrameState &frame) {
    if constexpr (Outer)
      return frame.prepared_outer;
    else
      return frame.prepared_inner;
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
 * @tparam BindingT The look's Binding.
 * @tparam TrackPath Whether the stage accumulates path length.
 * @pre The look's `surface_type` is SurfaceNoiseParams, so the frame carries a
 *      noise pointer and a loop offset.
 */
template <typename BindingT, bool TrackPath = false> struct SurfaceProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.surface.noise;
  }
  static const Vector &loop_offset(const FrameState &frame) {
    return frame.surface.loop_offset;
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
 * accessors a look's chosen source policy names are instantiated.
 */
template <typename BindingT> struct SourceProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static const auto &params(const FrameState &frame) {
    return frame.params.source;
  }
  static const PreparedSource &prepared(const FrameState &frame) {
    return frame.prepared_source;
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
 * @details Names all three fields the value families define; a look's material
 * stage instantiates only the accessors its transfer and coverage policies
 * call, so an IsoValueParams look never touches `edge_width` and vice versa.
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
 * @tparam BindingT The look's Binding.
 * @tparam HueV Hue-rotation source, translated to the pullback enum.
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
  static Pullback::Color::HueMode hue_mode(const FrameState &) {
    if constexpr (HueV == HueMode::NOISE)
      return Pullback::Color::HueMode::NOISE;
    if constexpr (HueV == HueMode::PATH_LENGTH)
      return Pullback::Color::HueMode::PATH_LENGTH;
    return Pullback::Color::HueMode::NONE;
  }
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

/** @brief Linear interpolation with exact endpoints, for a scalar field. */
inline float lerp(float from, float to, float progress) {
  return FixedPipeline::linear(from, to, progress);
}

/**
 * @brief Interpolates one parameter family across a preset transition.
 * @details One overload per family, each picking the interpolator its fields'
 * domains call for: linear for unconstrained scalars, geometric for positive
 * scales, shortest-arc for wrapping angles. The non-interpolable fields, namely
 * the Mobius coefficients and the palette mapping enum, snap to @p b once
 * progress reaches 1.
 * @param a Value at progress 0.
 * @param b Value at progress 1.
 * @param t Progress fraction.
 * @return The interpolated family.
 */
inline GridSourceParams interpolate(const GridSourceParams &a,
                                    const GridSourceParams &b, float t) {
  return {lerp(a.pattern_freq, b.pattern_freq, t),
          lerp(a.speed, b.speed, t),
          lerp(a.complexity, b.complexity, t),
          lerp(a.pattern_mix, b.pattern_mix, t),
          lerp(a.secondary_rate, b.secondary_rate, t),
          lerp(a.angle_rate, b.angle_rate, t)};
}

inline TwinWaveSourceParams interpolate(const TwinWaveSourceParams &a,
                                        const TwinWaveSourceParams &b,
                                        float t) {
  return {lerp(a.pattern_freq, b.pattern_freq, t), lerp(a.speed, b.speed, t),
          lerp(a.secondary_rate, b.secondary_rate, t),
          lerp(a.angle_rate, b.angle_rate, t)};
}

inline NoiseSourceParams interpolate(const NoiseSourceParams &a,
                                     const NoiseSourceParams &b, float t) {
  return {lerp(a.noise_scale, b.noise_scale, t),
          lerp(a.noise_contrast, b.noise_contrast, t),
          lerp(a.noise_time_rate, b.noise_time_rate, t)};
}

inline LatticeSourceParams interpolate(const LatticeSourceParams &a,
                                       const LatticeSourceParams &b, float t) {
  return {
      FixedPipeline::log_positive(a.lattice_cell_scale, b.lattice_cell_scale,
                                  t),
      lerp(a.lattice_shape_blend, b.lattice_shape_blend, t),
      FixedPipeline::log_positive(a.lattice_softness, b.lattice_softness, t),
      lerp(a.lattice_radius, b.lattice_radius, t)};
}

inline NoSurfaceParams interpolate(const NoSurfaceParams &,
                                   const NoSurfaceParams &, float) {
  return {};
}

inline SurfaceNoiseParams interpolate(const SurfaceNoiseParams &a,
                                      const SurfaceNoiseParams &b, float t) {
  return {lerp(a.scale, b.scale, t), lerp(a.strength, b.strength, t),
          lerp(a.speed, b.speed, t)};
}

inline ProjectionParams interpolate(const ProjectionParams &a,
                                    const ProjectionParams &b, float t) {
  return {lerp(a.pole_fade, b.pole_fade, t), lerp(a.spin_rate, b.spin_rate, t),
          lerp(a.wander, b.wander, t),
          lerp(a.camera_wander, b.camera_wander, t),
          FixedPipeline::shortest_periodic(a.central_meridian,
                                           b.central_meridian, t, TWO_PI_F)};
}

inline NoWarpParams interpolate(const NoWarpParams &a, const NoWarpParams &b,
                                float t) {
  return {lerp(a.speed, b.speed, t)};
}

inline MirrorParams interpolate(const MirrorParams &a, const MirrorParams &b,
                                float t) {
  return {lerp(a.speed, b.speed, t),
          FixedPipeline::shortest_periodic(a.rotation, b.rotation, t, TWO_PI_F),
          FixedPipeline::log_positive(a.cell_x, b.cell_x, t),
          FixedPipeline::log_positive(a.cell_y, b.cell_y, t),
          lerp(a.offset_x, b.offset_x, t),
          lerp(a.offset_y, b.offset_y, t)};
}

inline WaveShearParams interpolate(const WaveShearParams &a,
                                   const WaveShearParams &b, float t) {
  return {lerp(a.speed, b.speed, t), lerp(a.strength, b.strength, t),
          FixedPipeline::log_positive(a.frequency, b.frequency, t),
          FixedPipeline::shortest_periodic(a.field_angle, b.field_angle, t,
                                           TWO_PI_F),
          lerp(a.edge_width, b.edge_width, t)};
}

inline VectorNoiseParams interpolate(const VectorNoiseParams &a,
                                     const VectorNoiseParams &b, float t) {
  return {lerp(a.speed, b.speed, t), lerp(a.strength, b.strength, t),
          FixedPipeline::log_positive(a.scale, b.scale, t),
          FixedPipeline::shortest_periodic(a.vector_angle, b.vector_angle, t,
                                           TWO_PI_F),
          lerp(a.edge_width, b.edge_width, t)};
}

inline AffineParams interpolate(const AffineParams &a, const AffineParams &b,
                                float t) {
  return {lerp(a.speed, b.speed, t),
          lerp(a.rotation_rate, b.rotation_rate, t),
          lerp(a.translation_x, b.translation_x, t),
          lerp(a.translation_y, b.translation_y, t),
          FixedPipeline::log_positive(a.scale_x, b.scale_x, t),
          FixedPipeline::log_positive(a.scale_y, b.scale_y, t),
          lerp(a.shear, b.shear, t)};
}

inline PolarParams interpolate(const PolarParams &a, const PolarParams &b,
                               float t) {
  return {lerp(a.speed, b.speed, t),
          FixedPipeline::log_positive(a.radial_scale, b.radial_scale, t),
          FixedPipeline::shortest_periodic(a.radial_phase, b.radial_phase, t,
                                           TWO_PI_F),
          FixedPipeline::shortest_periodic(a.angular_phase, b.angular_phase, t,
                                           TWO_PI_F)};
}

inline NoLensParams interpolate(const NoLensParams &, const NoLensParams &,
                                float) {
  return {};
}

inline MobiusLensParams interpolate(const MobiusLensParams &a,
                                    const MobiusLensParams &b, float t) {
  MobiusLensParams value;
  value.mobius = t < 1.0f ? a.mobius : b.mobius;
  return value;
}

inline LinearValueParams interpolate(const LinearValueParams &,
                                     const LinearValueParams &, float) {
  return {};
}

inline EdgeValueParams interpolate(const EdgeValueParams &a,
                                   const EdgeValueParams &b, float t) {
  return {lerp(a.edge_width, b.edge_width, t)};
}

inline IsoValueParams interpolate(const IsoValueParams &a,
                                  const IsoValueParams &b, float t) {
  return {lerp(a.iso_level, b.iso_level, t),
          FixedPipeline::log_positive(a.iso_width, b.iso_width, t)};
}

inline ColorParams interpolate(const ColorParams &a, const ColorParams &b,
                               float t) {
  return {
      lerp(a.hue_shift_amount, b.hue_shift_amount, t),
      FixedPipeline::log_positive(a.hue_noise_scale, b.hue_noise_scale, t),
      lerp(a.hue_noise_speed, b.hue_noise_speed, t),
      lerp(a.palette_chroma, b.palette_chroma, t),
      FixedPipeline::log_positive(a.mapping_frequency, b.mapping_frequency, t),
      lerp(a.mapping_phase, b.mapping_phase, t),
      lerp(a.phase_oscillation_depth, b.phase_oscillation_depth, t),
      lerp(a.phase_oscillation_speed, b.phase_oscillation_speed, t),
      lerp(a.brightness_depth, b.brightness_depth, t),
      lerp(a.opacity_low, b.opacity_low, t),
      lerp(a.opacity_high, b.opacity_high, t),
      t < 1.0f ? a.palette_mapping : b.palette_mapping};
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
 * @brief Whether a scalar is finite and inside an inclusive range.
 * @return False for NaN, which fails both comparisons.
 */
inline bool finite_range(float value, float minimum, float maximum) {
  return std::isfinite(value) && value >= minimum && value <= maximum;
}

/**
 * @brief Whether every field of a parameter family is inside its authored
 *        range.
 * @details One overload per family, each mirroring the bounds
 * `Runtime::register_parameters` gives that family's sliders. This is the
 * admissibility test a restored snapshot must pass.
 * @param p Family to check.
 * @return True when every field is finite and in range.
 */
inline bool valid(const GridSourceParams &p) {
  return finite_range(p.pattern_freq, 0.1f, 20.0f) &&
         finite_range(p.speed, 0.0f, 0.5f) &&
         finite_range(p.complexity, 0.0f, 3.0f) &&
         finite_range(p.pattern_mix, 0.0f, 1.0f) &&
         finite_range(p.secondary_rate, 0.0f, 1.25f) &&
         finite_range(p.angle_rate, 0.0f, 0.05f);
}

inline bool valid(const TwinWaveSourceParams &p) {
  return finite_range(p.pattern_freq, 0.1f, 20.0f) &&
         finite_range(p.speed, 0.0f, 0.5f) &&
         finite_range(p.secondary_rate, 0.0f, 1.25f) &&
         finite_range(p.angle_rate, 0.0f, 0.05f);
}

inline bool valid(const NoiseSourceParams &p) {
  return finite_range(p.noise_scale, 1.0f / 64.0f, 64.0f) &&
         finite_range(p.noise_contrast, 0.0f, 8.0f) &&
         finite_range(p.noise_time_rate, -1.0f / 64.0f, 1.0f / 64.0f);
}

inline bool valid(const LatticeSourceParams &p) {
  return finite_range(p.lattice_cell_scale, 1.0f / 64.0f, 8.0f) &&
         finite_range(p.lattice_shape_blend, 0.0f, 1.0f) &&
         finite_range(p.lattice_softness, 1.0f / 1024.0f, 1.0f) &&
         finite_range(p.lattice_radius, 1.0f / 64.0f, 0.49f);
}

inline bool valid(const NoSurfaceParams &) { return true; }

inline bool valid(const SurfaceNoiseParams &p) {
  return finite_range(p.scale, 1.0f / 64.0f, 64.0f) &&
         finite_range(p.strength, -0.5f, 0.5f) &&
         finite_range(p.speed, -1.0f / 64.0f, 1.0f / 64.0f);
}

inline bool valid(const ProjectionParams &p) {
  return finite_range(p.pole_fade, 1.0f, 20.0f) &&
         finite_range(p.spin_rate, 0.0f, 0.05f) &&
         finite_range(p.wander, 0.0f, 1.0f) &&
         finite_range(p.camera_wander, 0.0f, 1.0f) &&
         finite_range(p.central_meridian, 0.0f, TWO_PI_F);
}

inline bool valid(const NoWarpParams &p) {
  return finite_range(p.speed, -0.02f, 0.02f);
}

inline bool valid(const MirrorParams &p) {
  return finite_range(p.speed, -0.02f, 0.02f) &&
         finite_range(p.rotation, 0.0f, TWO_PI_F) &&
         finite_range(p.cell_x, 1.0f / 64.0f, 8.0f) &&
         finite_range(p.cell_y, 1.0f / 64.0f, 8.0f) &&
         finite_range(p.offset_x, -8.0f, 8.0f) &&
         finite_range(p.offset_y, -8.0f, 8.0f);
}

inline bool valid(const WaveShearParams &p) {
  return finite_range(p.speed, -0.02f, 0.02f) &&
         finite_range(p.strength, -30.0f, 30.0f) &&
         finite_range(p.frequency, 0.01f, 32.0f) &&
         finite_range(p.field_angle, 0.0f, TWO_PI_F) &&
         finite_range(p.edge_width, 0.0f, 1.0f);
}

inline bool valid(const VectorNoiseParams &p) {
  return finite_range(p.speed, -0.02f, 0.02f) &&
         finite_range(p.strength, -30.0f, 30.0f) &&
         finite_range(p.scale, 1.0f / 64.0f, 64.0f) &&
         finite_range(p.vector_angle, 0.0f, TWO_PI_F) &&
         finite_range(p.edge_width, 0.0f, 1.0f);
}

inline bool valid(const AffineParams &p) {
  return finite_range(p.speed, -0.02f, 0.02f) &&
         finite_range(p.rotation_rate, -TWO_PI_F, TWO_PI_F) &&
         finite_range(p.translation_x, -4.0f, 4.0f) &&
         finite_range(p.translation_y, -4.0f, 4.0f) &&
         finite_range(p.scale_x, 1.0f / 64.0f, 64.0f) &&
         finite_range(p.scale_y, 1.0f / 64.0f, 64.0f) &&
         finite_range(p.shear, -4.0f, 4.0f);
}

inline bool valid(const PolarParams &p) {
  return finite_range(p.speed, -0.02f, 0.02f) &&
         finite_range(p.radial_scale, 1.0f / 64.0f, 64.0f) &&
         finite_range(p.radial_phase, -TWO_PI_F, TWO_PI_F) &&
         finite_range(p.angular_phase, -TWO_PI_F, TWO_PI_F);
}

inline bool valid(const NoLensParams &) { return true; }
inline bool valid(const MobiusLensParams &p) {
  const float values[] = {p.mobius.a.re, p.mobius.a.im, p.mobius.b.re,
                          p.mobius.b.im, p.mobius.c.re, p.mobius.c.im,
                          p.mobius.d.re, p.mobius.d.im};
  for (float value : values)
    if (!std::isfinite(value))
      return false;
  return true;
}
inline bool valid(const LinearValueParams &) { return true; }
inline bool valid(const EdgeValueParams &p) {
  return finite_range(p.edge_width, 0.0f, 1.0f);
}
inline bool valid(const IsoValueParams &p) {
  return finite_range(p.iso_level, 0.0f, 1.0f) &&
         finite_range(p.iso_width, 1.0f / 1024.0f, 1.0f);
}
inline bool valid(const ColorParams &p) {
  return finite_range(p.hue_shift_amount, -4.0f, 4.0f) &&
         finite_range(p.hue_noise_scale, 1.0f / 64.0f, 8.0f) &&
         finite_range(p.hue_noise_speed, -0.001f, 0.001f) &&
         finite_range(p.palette_chroma, 0.0f, 1.0f) &&
         finite_range(p.mapping_frequency, 1.0f, 32.0f) &&
         finite_range(p.mapping_phase, -1.0f, 1.0f) &&
         finite_range(p.phase_oscillation_depth, 0.0f, 1.0f) &&
         finite_range(p.phase_oscillation_speed, -0.01f, 0.01f) &&
         finite_range(p.brightness_depth, 0.0f, 1.0f) &&
         finite_range(p.opacity_low, 0.0f, 1.0f) &&
         finite_range(p.opacity_high, 0.0f, 1.0f) &&
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

/**
 * @brief Shared lifecycle for the fixed-look effects: parameter registration,
 * preset choreography, palette cycling, camera walks and noise clocks.
 * @details Curiously recurring: `Derived` supplies the look (its render
 * pipeline and preset bank) and inherits everything else. Required `Derived`
 * members are `PRESET_IDS`, `PARAMETER_SCHEMA_VERSION`, `PRESET_DWELL_FRAMES`,
 * `ANIMATED_PROJECTION`, `shade`, `initial_params` and `preset_params`, plus an
 * `OUTER_NOISE_SEED` / `SOURCE_NOISE_SEED` / `SURFACE_NOISE_SEED` for each
 * noise field the parameter set and `Has*Noise` flags request. Optional members,
 * detected by `requires` and defaulted when absent, are `ANIMATED_MOBIUS`,
 * `USES_CENTRAL_MERIDIAN` and an `after_fixed_init()` hook.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam Derived The effect class deriving from this runtime.
 * @tparam ParamsT The effect's `FixedLook::Params` specialization.
 * @tparam Harmony Palette harmony the generated palettes are drawn from.
 * @tparam HueV Hue-rotation source: none, noise field, or path length.
 * @tparam BrightnessV Brightness envelope applied by the color stage.
 * @tparam HasOuterNoise Whether the outer-camera stage owns a noise field.
 * @tparam HasSourceNoise Whether the source stage owns a noise field.
 */
template <int W, int H, typename Derived, typename ParamsT,
          PaletteHarmony Harmony, HueMode HueV,
          Pullback::Color::BrightnessEnvelope BrightnessV,
          bool HasOuterNoise = false, bool HasSourceNoise = false>
class Runtime : public Effect {
public:
  using Params = ParamsT;
  using Frame = FrameState<Params>;
  using PipelineBinding = Binding<Frame>;
  /** Frames an automatic preset transition spans. */
  static constexpr uint16_t TRANSITION_DURATION = 480;
  /** Whether the look displaces its surface, and so owns a surface noise
      field and needs a `Derived::SURFACE_NOISE_SEED`. */
  static constexpr bool HAS_SURFACE_NOISE =
      !std::is_same_v<typename ParamsT::surface_type, NoSurfaceParams>;

  /** @brief Constructs the effect at W x H with the POV column strobe on. */
  HS_COLD_MEMBER Runtime() : Effect(W, H, {.strobe = true}) {}

  /**
   * @brief Claims the runtime's persistent storage and registers the look.
   * @details Every allocation this makes comes from `persistent_arena`, so the
   * effect heap-allocates nothing after init. `Derived::after_fixed_init()`
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
    if constexpr (requires(Derived &effect) { effect.after_fixed_init(); })
      static_cast<Derived &>(*this).after_fixed_init();
  }

  /**
   * @brief Advances the frame clocks, resolves the frame state and shades.
   * @details The order is load-bearing: the runtime's clocks, camera walks and
   * palette all step before prepare_frame() snapshots them, so the whole scan
   * shades from one consistent frame. The transition evaluation is retired
   * after the scan, leaving the shaded frame and the transition's progress in
   * step.
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
      prepare_transition_value();
      advance_runtime();
      update_spatial_frames();
      update_palette_chroma();
      palette_cycler.step();
    }
    const Frame frame = prepare_frame();
    {
      HS_PROFILE(fx_shader_draw);
      Scan::Shader::draw_cached<W, H, 1>(canvas, [&frame](const Vector &view) {
        return Derived::shade(view, frame);
      });
    }
    finish_transition_evaluation();
  }

  /** @brief A parameter set tagged with the schema version that produced it. */
  struct ParameterSnapshot {
    /** `Derived::PARAMETER_SCHEMA_VERSION` at capture time. */
    uint32_t schema_version;
    Params params; /**< The captured parameter values. */
  };

  /**
   * @brief Captures the live parameters with the look's schema version.
   * @return A snapshot restore_parameters() accepts while the look's schema
   *         version is unchanged.
   */
  ParameterSnapshot serialize_parameters() const {
    return {Derived::PARAMETER_SCHEMA_VERSION, params};
  }

  /**
   * @brief Adopts a snapshot's parameters if it is admissible.
   * @details A look bumps `PARAMETER_SCHEMA_VERSION` whenever its `Params`
   * layout changes, so a snapshot taken under a different layout is rejected
   * rather than reinterpreted. On success any in-flight preset transition is
   * cancelled and the preset dwell restarts; on rejection nothing is touched.
   * @param snapshot Snapshot to adopt.
   * @return True when the snapshot's schema version matches and its parameters
   *         are valid, false otherwise.
   */
  bool restore_parameters(const ParameterSnapshot &snapshot) {
    if (snapshot.schema_version != Derived::PARAMETER_SCHEMA_VERSION ||
        !FixedLook::valid(snapshot.params))
      return false;
    transition.active = false;
    params = snapshot.params;
    palette_mapping = Pullback::Color::PaletteMappingWeights::single(
        params.color.palette_mapping);
    preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
    return true;
  }

protected:
  /** Descriptors the arena-backed parameter array holds; every slider a look
      registers has to fit. */
  static constexpr size_t PARAM_CAPACITY = 48;

  /** @brief One authored preset's parameter set. */
  struct Preset {
    Params params;
  };

  /** @brief An automatic preset transition and how far along it is. */
  struct Transition {
    Params from{}; /**< Parameters the transition departs from. */
    Params to{};   /**< Parameters it lands on. */
    Pullback::Color::PaletteMappingWeights mapping_from{};
    Pullback::Color::PaletteMappingWeights mapping_to{};
    /** Evaluations elapsed, counting to `duration`. */
    uint16_t evaluation = 0;
    uint16_t duration = TRANSITION_DURATION;
    bool active = false; /**< False when no transition is in flight. */
  };

  /**
   * @brief Starts a repeating circular Mobius warp on the lens parameters.
   * @details Pausable, so the warp freezes with the rest of the parameter
   * animations. Compiles to nothing for a look whose lens family carries no
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

  /**
   * @brief Overrides the dwell remaining before the next automatic transition.
   * @param frames Frames to hold; 0 lets the next frame start a transition.
   */
  HS_COLD_MEMBER void hold_initial_preset(uint16_t frames) {
    preset_dwell_remaining = frames;
  }

  /**
   * @brief Sets this frame's parameters from the in-flight transition's eased
   *        progress.
   * @details Under `Derived::ANIMATED_MOBIUS` the live lens coefficients are
   * carried across the interpolation, so the transition does not overwrite the
   * timeline animation driving them.
   */
  HS_COLD_MEMBER void prepare_transition_value() {
    if (!transition.active)
      return;
    MobiusParams animated_mobius;
    if constexpr (requires { Derived::ANIMATED_MOBIUS; })
      if constexpr (Derived::ANIMATED_MOBIUS)
        animated_mobius = params.lens.mobius;
    const FixedPipeline::EdgeProgress progress =
        FixedPipeline::edge_progress(transition.evaluation, transition.duration,
                                     FixedPipeline::Easing::EASE_IN_OUT_SIN);
    params =
        FixedLook::interpolate(transition.from, transition.to, progress.eased);
    if constexpr (requires { Derived::ANIMATED_MOBIUS; })
      if constexpr (Derived::ANIMATED_MOBIUS)
        params.lens.mobius = animated_mobius;
    palette_mapping = Pullback::Color::PaletteMappingWeights::lerp(
        transition.mapping_from, transition.mapping_to, progress.eased);
  }

  /// Advances the in-flight preset transition by one evaluation, latching the
  /// destination params once it reaches its duration.
  /// @details Deliberately not pause-gated: a started transition always runs to
  /// its authored endpoint so the look never rests between two presets.
  HS_COLD_MEMBER void finish_transition_evaluation() {
    if (!transition.active)
      return;
    if (transition.evaluation == transition.duration) {
      MobiusParams animated_mobius;
      if constexpr (requires { Derived::ANIMATED_MOBIUS; })
        if constexpr (Derived::ANIMATED_MOBIUS)
          animated_mobius = params.lens.mobius;
      params = transition.to;
      if constexpr (requires { Derived::ANIMATED_MOBIUS; })
        if constexpr (Derived::ANIMATED_MOBIUS)
          params.lens.mobius = animated_mobius;
      palette_mapping = transition.mapping_to;
      transition.active = false;
      preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
      return;
    }
    ++transition.evaluation;
  }

  /**
   * @brief Adopts a preset, crossfading only when choreography asked for it.
   * @details An AUTOMATIC change arms a TRANSITION_DURATION crossfade from the
   * live parameters; a manual or synchronized change snaps, since a user
   * driving the control expects the preset it names immediately.
   * @param change The requested preset move.
   * @return Always true: the runtime never vetoes a preset.
   */
  HS_COLD_MEMBER bool apply_preset(const PresetChange &change) override {
    const Params target = Derived::preset_params(change.to);
    if (change.origin != PresetChangeOrigin::AUTOMATIC) {
      transition.active = false;
      params = target;
      palette_mapping = Pullback::Color::PaletteMappingWeights::single(
          target.color.palette_mapping);
      preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
      return true;
    }
    transition = {params,
                  target,
                  palette_mapping,
                  Pullback::Color::PaletteMappingWeights::single(
                      target.color.palette_mapping),
                  0,
                  TRANSITION_DURATION,
                  true};
    return true;
  }

  /** Live parameters; the registered sliders write straight into these. */
  Params params = Derived::initial_params();
  Transition transition;

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
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "FixedLook::Runtime persistent footprint exceeds the default "
                "partition");

  static void configure_noise(FastNoiseLite &noise, int32_t seed) {
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise.SetSeed(seed);
    noise.SetFrequency(1.0f);
  }

  template <typename T> void register_source(T &source) {
    if constexpr (requires { source.pattern_freq; })
      register_animated_param("Pattern Freq", &source.pattern_freq, 0.1f,
                              20.0f);
    if constexpr (requires { source.speed; })
      register_animated_param("Speed", &source.speed, 0.0f, 0.5f);
    if constexpr (requires { source.complexity; })
      register_animated_param("Complexity", &source.complexity, 0.0f, 3.0f);
    if constexpr (requires { source.pattern_mix; })
      register_animated_param("Pattern Mix", &source.pattern_mix, 0.0f, 1.0f);
    if constexpr (requires { source.secondary_rate; })
      register_animated_param("Drift", &source.secondary_rate, 0.0f, 1.25f);
    if constexpr (requires { source.angle_rate; })
      register_animated_param("Source Angle Speed", &source.angle_rate, 0.0f,
                              0.05f);
    if constexpr (requires { source.noise_scale; })
      register_animated_param("Source Noise Scale", &source.noise_scale,
                              1.0f / 64.0f, 64.0f);
    if constexpr (requires { source.noise_contrast; })
      register_animated_param("Source Noise Contrast", &source.noise_contrast,
                              0.0f, 8.0f);
    if constexpr (requires { source.noise_time_rate; })
      register_animated_param("Source Noise Speed", &source.noise_time_rate,
                              -1.0f / 64.0f, 1.0f / 64.0f);
    if constexpr (requires { source.lattice_cell_scale; }) {
      register_animated_param("Lattice Cell Scale", &source.lattice_cell_scale,
                              1.0f / 64.0f, 8.0f);
      register_animated_param("Lattice Shape", &source.lattice_shape_blend,
                              0.0f, 1.0f);
      register_animated_param("Lattice Softness", &source.lattice_softness,
                              1.0f / 1024.0f, 1.0f);
      register_animated_param("Lattice Radius", &source.lattice_radius,
                              1.0f / 64.0f, 0.49f);
    }
  }

  template <typename T> void register_surface(T &surface) {
    if constexpr (!std::is_same_v<T, NoSurfaceParams>) {
      register_animated_param("Surface Noise Scale", &surface.scale,
                              1.0f / 64.0f, 64.0f);
      register_animated_param("Surface Noise Strength", &surface.strength,
                              -0.5f, 0.5f);
      register_animated_param("Surface Noise Speed", &surface.speed,
                              -1.0f / 64.0f, 1.0f / 64.0f);
    }
  }

  template <typename T> void register_warp(T &warp, const char *prefix) {
    if constexpr (!std::is_same_v<T, NoWarpParams>) {
      register_animated_param(prefix, &warp.speed, -0.02f, 0.02f);
      if constexpr (requires { warp.strength; })
        register_animated_param("Warp Strength", &warp.strength, -30.0f, 30.0f);
      if constexpr (requires { warp.frequency; })
        register_animated_param("Warp Frequency", &warp.frequency, 0.01f,
                                32.0f);
      if constexpr (requires { warp.field_angle; })
        register_animated_param("Warp Field Angle", &warp.field_angle, 0.0f,
                                TWO_PI_F);
      if constexpr (requires { warp.scale; })
        register_animated_param("Warp Scale", &warp.scale, 1.0f / 64.0f, 64.0f);
      if constexpr (requires { warp.vector_angle; })
        register_animated_param("Warp Vector Angle", &warp.vector_angle, 0.0f,
                                TWO_PI_F);
      if constexpr (std::is_same_v<T, MirrorParams>) {
        register_animated_param("Mirror Rotation", &warp.rotation, 0.0f,
                                TWO_PI_F);
        register_animated_param("Mirror Cell X", &warp.cell_x, 1.0f / 64.0f,
                                8.0f);
        register_animated_param("Mirror Cell Y", &warp.cell_y, 1.0f / 64.0f,
                                8.0f);
        register_animated_param("Mirror Offset X", &warp.offset_x, -8.0f, 8.0f);
        register_animated_param("Mirror Offset Y", &warp.offset_y, -8.0f, 8.0f);
      }
      if constexpr (std::is_same_v<T, AffineParams>) {
        register_animated_param("Affine Rotation Rate", &warp.rotation_rate,
                                -TWO_PI_F, TWO_PI_F);
        register_animated_param("Affine Translation X", &warp.translation_x,
                                -4.0f, 4.0f);
        register_animated_param("Affine Translation Y", &warp.translation_y,
                                -4.0f, 4.0f);
        register_animated_param("Affine Scale X", &warp.scale_x, 1.0f / 64.0f,
                                64.0f);
        register_animated_param("Affine Scale Y", &warp.scale_y, 1.0f / 64.0f,
                                64.0f);
        register_animated_param("Affine Shear", &warp.shear, -4.0f, 4.0f);
      }
      if constexpr (std::is_same_v<T, PolarParams>) {
        register_animated_param("Polar Radial Scale", &warp.radial_scale,
                                1.0f / 64.0f, 64.0f);
        register_animated_param("Polar Radial Phase", &warp.radial_phase,
                                -TWO_PI_F, TWO_PI_F);
        register_animated_param("Polar Angular Phase", &warp.angular_phase,
                                -TWO_PI_F, TWO_PI_F);
      }
    }
  }

  void register_parameters() {
    register_source(params.source);
    register_animated_param("Pole Fade", &params.projection.pole_fade, 1.0f,
                            20.0f);
    if constexpr (Derived::ANIMATED_PROJECTION) {
      register_animated_param("Projection Spin Speed",
                              &params.projection.spin_rate, 0.0f, 0.05f);
      register_animated_param("Projection Wander", &params.projection.wander,
                              0.0f, 1.0f);
    }
    register_animated_param("Camera Wander", &params.projection.camera_wander,
                            0.0f, 1.0f);
    if constexpr (requires { Derived::USES_CENTRAL_MERIDIAN; })
      if constexpr (Derived::USES_CENTRAL_MERIDIAN)
        register_animated_param("Central Meridian",
                                &params.projection.central_meridian, 0.0f,
                                TWO_PI_F);
    register_surface(params.surface);
    register_warp(params.outer_warp, "Planar Warp 1 Speed");
    register_warp(params.inner_warp, "Planar Warp 2 Speed");
    if constexpr (requires { params.value.edge_width; })
      register_animated_param("Edge Width", &params.value.edge_width, 0.0f,
                              1.0f);
    if constexpr (requires { params.value.iso_level; }) {
      register_animated_param("Iso Level", &params.value.iso_level, 0.0f, 1.0f);
      register_animated_param("Iso Width", &params.value.iso_width,
                              1.0f / 1024.0f, 1.0f);
    }
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

  /// Retires the preset dwell and starts the next automatic preset transition.
  /// @details Pause suppresses preset selection, so no new transition begins
  /// while paused.
  HS_COLD_MEMBER void begin_automatic_transition() {
    if constexpr (Derived::PRESET_IDS.size() == 1)
      return;
    if (anims_paused || transition.active)
      return;
    if (preset_dwell_remaining > 0 && --preset_dwell_remaining > 0)
      return;
    if (advancePreset()) {
#ifdef HS_PROFILE_ENABLE
      hs::log("Preset: %u/%u", static_cast<unsigned>(getPresetIndex() + 1),
              static_cast<unsigned>(getPresetCount()));
#endif
    }
  }

  /**
   * @brief Steps every phase clock the look's parameter families define.
   * @details Each clock is compiled in only when its field exists, so a look
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
      projection_wander =
          (FixedPipeline::scaled_rotation_delta(projection_delta.normalized(),
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
    outer_wander =
        (FixedPipeline::scaled_rotation_delta(outer_delta.normalized(),
                                              params.projection.camera_wander) *
         outer_wander)
            .normalized();
    outer_conjugate = outer_wander.conjugate();
  }

  /**
   * @brief Resolves one warp slot's per-frame rotation and transform.
   * @details Only the affine slot rotates with the frame, and only in the
   * outer position; every other family takes its rotation straight from a
   * parameter.
   * @param warp The slot's parameters.
   * @param phase The slot's phase clock.
   * @param outer True for the first warp slot.
   * @return The slot's PreparedWarp, with the union member its family uses.
   */
  template <typename WarpT>
  HS_COLD_MEMBER PreparedWarp prepare_warp(const WarpT &warp, float phase,
                                           bool outer) const {
    PreparedWarp prepared{};
    float rotation = 0.0f;
    if constexpr (std::is_same_v<WarpT, MirrorParams>) {
      rotation = warp.rotation;
      prepared.transform.mirror = {
          wrap_t(warp.offset_x / warp.cell_x + phase) * warp.cell_x,
          wrap_t(warp.offset_y / warp.cell_y) * warp.cell_y};
    } else if constexpr (std::is_same_v<WarpT, WaveShearParams>) {
      rotation = warp.field_angle;
    } else if constexpr (std::is_same_v<WarpT, VectorNoiseParams>) {
      rotation = warp.vector_angle;
      const float angle = TWO_PI_F * wrap_t(phase);
      prepared.transform.noise_loop = {NOISE_LOOP_RADIUS * sinf(angle) *
                                           0.7071067811865475f,
                                       NOISE_LOOP_RADIUS * cosf(angle)};
    } else if constexpr (std::is_same_v<WarpT, AffineParams>) {
      static_assert(
          std::is_same_v<typename Params::source_type, LatticeSourceParams>,
          "the affine warp stage translates in lattice cells and requires a "
          "LatticeSourceParams source");
      const float cycle = TWO_PI_F * wrap_t(phase);
      const float cycle_cos = cosf(cycle);
      const float period = 1.0f / params.source.lattice_cell_scale;
      rotation = outer ? outer_rotation : 0.0f;
      prepared.transform.affine = {wrap_t(phase) * warp.translation_x * period,
                                   wrap_t(phase) * warp.translation_y * period,
                                   powf(warp.scale_x, cycle_cos),
                                   powf(warp.scale_y, cycle_cos),
                                   warp.shear * cycle_cos};
    }
    prepared.rotation_cos = cosf(rotation);
    prepared.rotation_sin = sinf(rotation);
    return prepared;
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
    PreparedSurface<typename Params::surface_type> surface{};
    if constexpr (HAS_SURFACE_NOISE) {
      const float angle = TWO_PI_F * wrap_t(surface_phase);
      surface = {&state->surface.noise,
                 Vector(NOISE_LOOP_RADIUS * cosf(angle),
                        NOISE_LOOP_RADIUS * sinf(angle), 0.0f)};
    }
    return {Derived::ANIMATED_PROJECTION ? projection_conjugate : Quaternion(),
            outer_conjugate,
            {source_primary, source_secondary, source_angle,
             fast_cosf(source_angle), fast_sinf(source_angle)},
            prepare_warp(params.outer_warp, outer_phase, true),
            prepare_warp(params.inner_warp, inner_phase, false),
            outer_noise,
            source_noise,
            &palette_cycler.palette(),
            state->hue_rotation_lut.data(),
            state->hue_noise_lut.data(),
            params,
            palette_mapping,
            outer_phase,
            inner_phase,
            source_noise_time,
            palette_oscillation_phase,
            surface};
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
   * palette is the one the look's authored hue names.
   * @param context The Runtime instance, as registered with the cycler.
   * @param sequence Zero-based index of the palette being generated.
   * @param out Receives the recipe to bake.
   */
  static void next_palette(void *context, uint32_t sequence,
                           GenerativePalette &out) {
    Runtime &effect = *static_cast<Runtime *>(context);
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
  uint16_t preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
  Timeline timeline;
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

} // namespace FixedLook
