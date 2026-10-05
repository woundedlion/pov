/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

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
