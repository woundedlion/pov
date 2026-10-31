/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/** @file composed_providers.h
 * @brief Providers binding pullback stages to a composed effect's
 * frame state and parameter families.
 */

/** @brief Pipeline order of the surface stage relative to the lens stage. */
enum class SurfacePlacement : uint8_t { BEFORE_LENS, AFTER_LENS };

// The stage vocabulary a composed effect names, re-exported from the
// per-stage headers.
using Color::ColorParams;
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
 * @details Its pointers alias the runtime's persistent state and the palette
 * cycler's current bake, so a frame stays valid only until the next
 * draw_frame().
 */
template <typename ParamsT> struct FrameState {
  /** Conjugate of the projection orientation; identity unless the effect sets
      `ANIMATED_PROJECTION`. */
  math::Quaternion projection_conjugate;
  /** Conjugate of the outer camera orientation. */
  math::Quaternion outer_conjugate;
  const BakedPalette *palette; /**< The cycler's current bake. */
  /** Current palette hue table. */
  const Pixel *hue_rotation_lut;
  ParamsT params; /**< The frame's parameter values, already interpolated if a
                       preset transition is in flight. */
  /** Palette mapping weights, blended across a preset transition. */
  Pullback::Color::PaletteMappingWeights palette_mapping;
  ComposedDetail::ResourceStorage<ComposedDetail::ResourceFrame,
                                  typename ParamsT::ResourceTypes>
      resources; ///< Per-resource clocks and noise views, keyed by resource.
  float palette_oscillation_phase; /**< Phase of the mapping wobble. */
  const PaletteLabLut *palette_lab_lut =
      nullptr; ///< Borrowed interval base cache.
};

/** @brief Column-sampled context with retained color resources. */
template <typename ParamsT> struct FrameIntervalState : FrameState<ParamsT> {
  Color::TemporalPaletteResources temporal; ///< Previous color resources.
};

/**
 * @brief Ties a pullback pipeline to one effect's frame state.
 * @details Composed effects render uninstrumented, so the binding pins
 * Pullback::NoInstrumentation.
 * @tparam FrameT The effect's FrameState specialization.
 */
template <typename FrameT> struct Binding {
  using FrameState = FrameT;                           ///< Frame type.
  using Instrumentation = Pullback::NoInstrumentation; ///< Uninstrumented.
};

/** @brief Supplies the camera orientation to Pullback::Stage::Rotate.
 * @tparam BindingT Frame binding of the effect. */
template <typename BindingT> struct OuterCameraProvider {
  using Binding = BindingT;                        ///< The effect's Binding.
  using FrameState = typename Binding::FrameState; ///< The effect's frame.
  /** @brief Conjugate of the outer camera orientation.
   *  @param frame Current frame.
   *  @return `frame.outer_conjugate`. */
  __attribute__((always_inline)) static const math::Quaternion &
  conjugate(const FrameState &frame) {
    return frame.outer_conjugate;
  }
};

/**
 * @brief Supplies the projection frame and its parameters to the
 *        Pullback::Projection policies.
 * @tparam BindingT Frame binding of the effect.
 */
template <typename BindingT> struct ProjectionProvider {
  using Binding = BindingT;                        ///< The effect's Binding.
  using FrameState = typename Binding::FrameState; ///< The effect's frame.
  /** @brief Conjugate of the projection orientation.
   *  @param frame Current frame.
   *  @return `frame.projection_conjugate`. */
  __attribute__((always_inline)) static const math::Quaternion &
  conjugate(const FrameState &frame) {
    return frame.projection_conjugate;
  }
  /** @brief Projection singularity fade width.
   *  @param frame Current frame.
   *  @return The projection family's `singularity_fade`. */
  __attribute__((always_inline)) static float
  singularity_fade(const FrameState &frame) {
    return frame.params.template get<"projection">().singularity_fade;
  }
  /** @brief Projection central meridian.
   *  @param frame Current frame.
   *  @return The projection family's `central_meridian`. */
  __attribute__((always_inline)) static float
  central_meridian(const FrameState &frame) {
    return frame.params.template get<"projection">().central_meridian;
  }
};

/** @brief Supplies the Mobius coefficients to Pullback::Lens::Mobius.
 * @tparam BindingT Frame binding of the effect.
 * @tparam Key Lens resource key. */
template <typename BindingT, ResourceKey Key = "lens"> struct LensProvider {
  using Binding = BindingT;                        ///< The effect's Binding.
  using FrameState = typename Binding::FrameState; ///< The effect's frame.
  /** @brief The lens's Mobius coefficients.
   *  @param frame Current frame.
   *  @return The `Key` family's `mobius`. */
  __attribute__((always_inline)) static const math::MobiusParams &
  params(const FrameState &frame) {
    return frame.params.template get<Key>().mobius;
  }
};

/** @brief Supplies parameters, clock and noise of one planar warp instance.
 * @tparam BindingT Frame binding of the effect.
 * @tparam Key Warp resource key.
 * @tparam Family Warp parameter family.
 * @tparam TrackPath Whether the stage accumulates path length.
 * @tparam SourceKey Source resource supplying affine lattice cell scale. */
template <typename BindingT, ResourceKey Key, typename Family,
          bool TrackPath = false, ResourceKey SourceKey = "source">
struct WarpProvider {
  using Binding = BindingT;                        ///< The effect's Binding.
  using FrameState = typename Binding::FrameState; ///< The effect's frame.
  /** @brief The warp instance's parameters.
   *  @param frame Current frame.
   *  @return The `Key` family. */
  __attribute__((always_inline)) static const auto &
  params(const FrameState &frame) {
    return frame.params.template get<Key>();
  }
  /** @brief Per-frame prepared warp state.
   *  @param frame Current frame.
   *  @return `Pullback::Warp::prepare` of the params and clocks, or
   *  `Pullback::NoPrepared` for a polar warp. */
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
  /** @brief The warp's phase clock.
   *  @param frame Current frame.
   *  @return Phase modulo 1; column-time samples may lie outside [0, 1). */
  __attribute__((always_inline)) static float phase(const FrameState &frame) {
    return frame.resources.template get<Key>().phase;
  }
  /** @brief The warp's noise generator.
   *  @param frame Current frame.
   *  @return The runtime-owned generator. */
  __attribute__((always_inline)) static const FastNoiseLite &
  noise(const FrameState &frame) {
    return *frame.resources.template get<Key>().noise;
  }
  /** @brief Whether the stage accumulates path length.
   *  @return `TrackPath`. */
  __attribute__((always_inline)) static bool
  path_length_required(const FrameState &) {
    return TrackPath;
  }
};

/**
 * @brief Supplies the displacement field to the Pullback::Surface policies.
 * @tparam BindingT The effect's Binding.
 * @tparam Family Displacement parameter family.
 * @tparam Key Surface resource key.
 * @tparam TrackPath Whether the stage accumulates path length.
 * @pre The resource family is a displacement family.
 */
template <typename BindingT, typename Family, bool TrackPath = false,
          ResourceKey Key = "surface">
struct SurfaceProvider {
  using Binding = BindingT;                        ///< The effect's Binding.
  using FrameState = typename Binding::FrameState; ///< The effect's frame.
  /** @brief The surface's noise generator.
   *  @param frame Current frame.
   *  @return The runtime-owned generator. */
  __attribute__((always_inline)) static const FastNoiseLite &
  noise(const FrameState &frame) {
    return *frame.resources.template get<Key>().noise;
  }
  /** @brief Per-frame prepared surface state.
   *  @param frame Current frame.
   *  @return `Pullback::Surface::prepare_direct` for a directed family,
   *  else `Pullback::Surface::prepare`. */
  __attribute__((always_inline)) static auto prepare(const FrameState &frame) {
    if constexpr (std::is_same_v<Family, SphereAnimation::BurstRippleParams>) {
      const auto &resource = frame.resources.template get<Key>();
      return SphereAnimation::prepare_ripples(params(frame), resource);
    } else if constexpr (requires {
                           frame.params.template get<Key>().direction;
                         })
      return Pullback::Surface::prepare_direct(
          frame.resources.template get<Key>().phase,
          frame.params.template get<Key>().direction);
    else
      return Pullback::Surface::prepare(
          frame.resources.template get<Key>().phase);
  }
  /** @brief Current-owned ripple curves selected for this column.
   * @param frame Current column context.
   * @return Transfer curves, or null for canonical evaluation. */
  __attribute__((always_inline)) static const SphereAnimation::RadialTables *
  radial_tables(const FrameState &frame) {
    const auto &resource = frame.resources.template get<Key>();
    return resource.radial_fraction >= 0.0f ? resource.radial_tables : nullptr;
  }
  /** @brief Ripple transfer weight of the current column.
   * @param frame Current column context.
   * @return Current endpoint weight. */
  __attribute__((always_inline)) static float
  radial_fraction(const FrameState &frame) {
    return frame.resources.template get<Key>().radial_fraction;
  }
  /** @brief The surface instance's parameters.
   *  @param frame Current frame.
   *  @return The `Key` family. */
  __attribute__((always_inline)) static const auto &
  params(const FrameState &frame) {
    return frame.params.template get<Key>();
  }
  /** @brief The surface's phase.
   *  @param frame Current frame.
   *  @return Ripple: the cycle in [0, 1); others: phase modulo 1,
   *  with column-time samples possibly outside [0, 1). */
  __attribute__((always_inline)) static float phase(const FrameState &frame) {
    using SurfaceParams = Family;
    if constexpr (std::is_same_v<SurfaceParams, PeriodicRippleParams>)
      return Pullback::Surface::ripple_cycle(
          frame.resources.template get<Key>().phase,
          frame.params.template get<Key>());
    else
      return frame.resources.template get<Key>().phase;
  }
  /** @brief Displacement field spatial scale.
   *  @param frame Current frame.
   *  @return The family's `scale`. */
  __attribute__((always_inline)) static float scale(const FrameState &frame) {
    return frame.params.template get<Key>().scale;
  }
  /** @brief Displacement strength.
   *  @param frame Current frame.
   *  @return The family's `strength`. */
  __attribute__((always_inline)) static float
  strength(const FrameState &frame) {
    return frame.params.template get<Key>().strength;
  }
  /** @brief Whether the stage accumulates path length.
   *  @return `TrackPath`. */
  __attribute__((always_inline)) static bool
  path_length_required(const FrameState &) {
    return TrackPath;
  }
};

namespace Surface {
/** @brief Sequential spherical ripple bursts.
 * @tparam State Surface provider. */
template <typename State> struct BurstRipple : ApproximationDefaults {
  using Binding = typename State::Binding;           ///< Frame binding.
  using FrameState = typename State::FrameState;     ///< Current frame.
  using Prepared = SphereAnimation::PreparedRipples; ///< Active waves.
  /** @brief Prepares the active wave kernels.
   * @param frame Current frame.
   * @return Active kernels in launch order. */
  static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }
  /** @brief Applies every active wave in launch order.
   * @param input Unit sphere direction.
   * @param frame Current frame.
   * @param prepared Active wave kernels.
   * @return Displaced direction and accumulated path length. */
  static SurfaceResult apply(const math::Vector &input, const FrameState &frame,
                             const Prepared &prepared) {
    if constexpr (requires {
                    State::radial_tables(frame);
                    State::radial_fraction(frame);
                  })
      return SphereAnimation::radial_ripples(
          input, prepared, State::path_length_required(frame),
          State::radial_tables(frame), State::radial_fraction(frame));
    else
      return SphereAnimation::grouped_ripples(
          input, prepared, State::path_length_required(frame));
  }
};
} // namespace Surface

/**
 * @brief Supplies the pattern and noise state to the Pullback::Source policies.
 * @details The pattern accessors read the advanced clocks. noise() and
 * noise_time() read runtime resources; noise_scale() and noise_contrast() read
 * the NoiseSourceParams fields.
 * @tparam BindingT Frame binding of the effect.
 * @tparam Family Source parameter family.
 * @tparam Key Source resource key.
 */
template <typename BindingT, typename Family, ResourceKey Key = "source">
struct SourceProvider {
  using Binding = BindingT;                        ///< The effect's Binding.
  using FrameState = typename Binding::FrameState; ///< The effect's frame.
  /** @brief The source instance's parameters.
   *  @param frame Current frame.
   *  @return The `Key` family. */
  __attribute__((always_inline)) static const auto &
  params(const FrameState &frame) {
    return frame.params.template get<Key>();
  }
  /** @brief Prepared pattern state from the source clocks.
   *  @param frame Current frame.
   *  @return `Pattern::prepare` of the primary, secondary and
   *  angle clocks. */
  __attribute__((always_inline)) static Pattern::PreparedSource
  prepare(const FrameState &frame) {
    return Pattern::prepare(frame.resources.template get<Key>().primary,
                            frame.resources.template get<Key>().secondary,
                            frame.resources.template get<Key>().angle);
  }
  /** @brief The source's noise generator.
   *  @param frame Current frame.
   *  @return The runtime-owned generator. */
  __attribute__((always_inline)) static const FastNoiseLite &
  noise(const FrameState &frame) {
    return *frame.resources.template get<Key>().noise;
  }
  /** @brief Noise spatial scale.
   *  @param frame Current frame.
   *  @return The family's `noise_scale`. */
  __attribute__((always_inline)) static float
  noise_scale(const FrameState &frame) {
    return frame.params.template get<Key>().noise_scale;
  }
  /** @brief Noise time coordinate.
   *  @param frame Current frame.
   *  @return The source clock's `noise_time` modulo 1; column-time samples
   *  may lie outside [0, 1). */
  __attribute__((always_inline)) static float
  noise_time(const FrameState &frame) {
    return frame.resources.template get<Key>().noise_time;
  }
  /** @brief Noise contrast.
   *  @param frame Current frame.
   *  @return The family's `noise_contrast`. */
  __attribute__((always_inline)) static float
  noise_contrast(const FrameState &frame) {
    return frame.params.template get<Key>().noise_contrast;
  }
};

/**
 * @brief Supplies the value-family fields to the Pullback::Transfer and
 *        Pullback::ValueCoverage and ProjectionCoverage::EdgeFade policies.
 * @tparam BindingT Frame binding of the effect.
 * @tparam Family Value parameter family.
 * @tparam Key Value resource key.
 */
template <typename BindingT, typename Family, ResourceKey Key = "value">
struct ValueProvider {
  using Binding = BindingT;                        ///< The effect's Binding.
  using FrameState = typename Binding::FrameState; ///< The effect's frame.
  /** @brief Iso-contour level.
   *  @param frame Current frame.
   *  @return The family's `iso_level`. */
  __attribute__((always_inline)) static float
  iso_level(const FrameState &frame) {
    static_assert(
        requires { frame.params.template get<Key>().iso_level; },
        "value provider requires iso_level");
    return frame.params.template get<Key>().iso_level;
  }
  /** @brief Iso-contour width.
   *  @param frame Current frame.
   *  @return The family's `iso_width`. */
  __attribute__((always_inline)) static float
  iso_width(const FrameState &frame) {
    static_assert(
        requires { frame.params.template get<Key>().iso_width; },
        "value provider requires iso_width");
    return frame.params.template get<Key>().iso_width;
  }
  /** @brief Projection edge-fade width.
   *  @param frame Current frame.
   *  @return The family's `edge_width`. */
  __attribute__((always_inline)) static float
  edge_width(const FrameState &frame) {
    static_assert(
        requires { frame.params.template get<Key>().edge_width; },
        "value provider requires edge_width");
    return frame.params.template get<Key>().edge_width;
  }
  /** @brief Value-cutout threshold.
   *  @param frame Current frame.
   *  @return The family's `cutout_threshold`. */
  __attribute__((always_inline)) static float
  cutout_threshold(const FrameState &frame) {
    static_assert(
        requires { frame.params.template get<Key>().cutout_threshold; },
        "value provider requires cutout_threshold");
    return frame.params.template get<Key>().cutout_threshold;
  }
  /** @brief Value-cutout softness.
   *  @param frame Current frame.
   *  @return The family's `cutout_softness`. */
  __attribute__((always_inline)) static float
  cutout_softness(const FrameState &frame) {
    static_assert(
        requires { frame.params.template get<Key>().cutout_softness; },
        "value provider requires cutout_softness");
    return frame.params.template get<Key>().cutout_softness;
  }
};

/**
 * @brief Supplies the palette, mapping and hue state to
 *        Pullback::Color::GeneratedPalette.
 * @details Source resources are immutable endpoint views.
 * @tparam BindingT The effect's Binding.
 * @tparam BrightnessV Brightness envelope reported to the color stage.
 * @tparam LabCacheV Whether the provider exposes an interval base cache.
 */
template <typename BindingT, Pullback::Color::BrightnessEnvelope BrightnessV,
          bool LabCacheV = false>
struct ColorProvider {
  using InteriorFractionOwner = ColorProvider; ///< Mirrors the frame fraction.
  using BrightnessOwner = ColorProvider; ///< Fixed-envelope contract owner.
  static constexpr auto BRIGHTNESS_ENVELOPE =
      BrightnessV;                                 ///< Fixed gain curve.
  using Binding = BindingT;                        ///< The effect's Binding.
  using FrameState = typename Binding::FrameState; ///< The effect's frame.
  /** @brief Borrows the interval's palette base cache.
   * @param frame Column frame context.
   * @return Persistent cache, or null for canonical draws. */
  static const PaletteLabLut *palette_lab_lut(const FrameState &frame)
    requires LabCacheV
  {
    return frame.palette_lab_lut;
  }
  /** @brief Previous color resources and the column's interpolation fraction.
   *  @param frame Current column context.
   *  @return The retained resources. */
  static Color::TemporalPaletteResources
  temporal_resources(const FrameState &frame)
    requires requires { frame.temporal; }
  {
    return frame.temporal;
  }
  /** @brief Palette mapping weights.
   *  @param frame Current frame.
   *  @return `frame.palette_mapping`. */
  __attribute__((always_inline)) static Pullback::Color::PaletteMappingWeights
  mapping_weights(const FrameState &frame) {
    return frame.palette_mapping;
  }
  /** @brief Palette mapping frequency.
   *  @param frame Current frame.
   *  @return The color family's `mapping_frequency`. */
  __attribute__((always_inline)) static float
  mapping_frequency(const FrameState &frame) {
    return frame.params.template get<"color">().mapping_frequency;
  }
  /** @brief Palette mapping phase.
   *  @param frame Current frame.
   *  @return The color family's `mapping_phase`. */
  __attribute__((always_inline)) static float
  mapping_phase(const FrameState &frame) {
    return frame.params.template get<"color">().mapping_phase;
  }
  /** @brief Mapping-phase oscillation depth.
   *  @param frame Current frame.
   *  @return The color family's `phase_oscillation_depth`. */
  __attribute__((always_inline)) static float
  oscillation_depth(const FrameState &frame) {
    return frame.params.template get<"color">().phase_oscillation_depth;
  }
  /** @brief Mapping-phase oscillation phase.
   *  @param frame Current frame.
   *  @return `frame.palette_oscillation_phase`. */
  __attribute__((always_inline)) static float
  oscillation_phase(const FrameState &frame) {
    return frame.palette_oscillation_phase;
  }
  /** @brief The cycler's current baked palette.
   *  @param frame Current frame.
   *  @return `*frame.palette`. */
  __attribute__((always_inline)) static const BakedPalette &
  palette(const FrameState &frame) {
    return *frame.palette;
  }
  /** @brief Hue-rotation LUT view.
   *  @param frame Current frame.
   *  @return The frame's palette hue table. */
  __attribute__((always_inline)) static Pullback::Color::HueRotationLutView
  hue_rotation(const FrameState &frame) {
    return {frame.hue_rotation_lut, frame.hue_rotation_lut != nullptr};
  }
  /** @brief Brightness envelope.
   *  @return `BrightnessV`. */
  __attribute__((always_inline)) static Pullback::Color::BrightnessEnvelope
  brightness_envelope(const FrameState &) {
    return BrightnessV;
  }
  /** @brief Brightness envelope floor.
   *  @param frame Current frame.
   *  @return The color family's `brightness_bottom`. */
  __attribute__((always_inline)) static float
  brightness_bottom(const FrameState &frame) {
    return frame.params.template get<"color">().brightness_bottom;
  }
  /** @brief Brightness envelope ceiling.
   *  @param frame Current frame.
   *  @return The color family's `brightness_top`. */
  __attribute__((always_inline)) static float
  brightness_top(const FrameState &frame) {
    return frame.params.template get<"color">().brightness_top;
  }
  /** @brief Opacity at value 0.
   *  @param frame Current frame.
   *  @return The color family's `opacity_low`. */
  __attribute__((always_inline)) static float
  opacity_low(const FrameState &frame) {
    return frame.params.template get<"color">().opacity_low;
  }
  /** @brief Opacity at value 1.
   *  @param frame Current frame.
   *  @return The color family's `opacity_high`. */
  __attribute__((always_inline)) static float
  opacity_high(const FrameState &frame) {
    return frame.params.template get<"color">().opacity_high;
  }
};

/** @brief Keyed source provider for an authored color modulation.
 * @tparam BindingT Pipeline binding.
 * @tparam Key Modulation instance key.
 * @tparam SourceV Source capability.
 * @tparam EvaluationV Source storage strategy.
 * @tparam PhaseV Whether eligible sources share column-phase evaluation. */
template <typename BindingT, ResourceKey Key, Pattern::SourceKind SourceV,
          Pattern::Evaluation EvaluationV = Pattern::Evaluation::BAKED_CUBE,
          bool PhaseV = false>
struct ModulationProvider {
  using Binding = BindingT;                        ///< Pipeline binding.
  using FrameState = typename Binding::FrameState; ///< Frame context.
  using Prepared = std::conditional_t<PhaseV, Color::PhaseModulationState,
                                      Color::ModulationState>; ///< Layer state.
  /** @brief Prepares endpoint views and mechanic controls.
   * @param frame Frame context.
   * @return Immutable layer state. */
  static Prepared prepare(const FrameState &frame) {
    const auto &config = frame.params.template get<Key>().modulation;
    const auto &resource = frame.resources.template get<Key>();
    const Color::ModulationState BASE{
        resource.modulation, resource.previous_modulation, config.source,
        config.mechanic, config.enabled};
    if constexpr (!PhaseV)
      return BASE;
    else {
      Color::PhaseModulationState result;
      static_cast<Color::ModulationState &>(result) = BASE;
      if constexpr (requires { frame.temporal.phase_compatible; }) {
        if (!frame.temporal.phase_compatible)
          return result;
        if (config.source != SourceV || config.evaluation != EvaluationV)
          return result;
        if (config.mechanic.kind != ColorModulation::Mechanic::HUE_ROTATION &&
            config.mechanic.kind !=
                ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT)
          return result;
        if constexpr (SourceV == Pattern::SourceKind::DISPLACEMENT)
          result.ready = true;
        else if constexpr (SourceV == Pattern::SourceKind::SPHERE_SMOOTH_GRID &&
                           EvaluationV == Pattern::Evaluation::DIRECT) {
          if (!resource.modulation || !resource.previous_modulation)
            return result;
          if (!Modulation::phase_compatible(
                  resource.previous_modulation->config,
                  resource.modulation->config))
            return result;
          const auto &a = std::get<3>(resource.previous_modulation->source);
          const auto &b = std::get<3>(resource.modulation->source);
          result.grid = b;
          result.grid.a.phase_turns =
              Spacetime::ClockInterval(a.a.phase_turns, b.a.phase_turns, 1.0f)
                  .sample(frame.temporal.fraction);
          result.grid.b.phase_turns =
              Spacetime::ClockInterval(a.b.phase_turns, b.b.phase_turns, 1.0f)
                  .sample(frame.temporal.fraction);
          result.ready = true;
        }
      }
      return result;
    }
  }
};
