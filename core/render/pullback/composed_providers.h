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
using Lens::MobiusGeneratorParams;
using Lens::TwistLensParams;
using SphereAnimation::RotationParams;
using Projection::ProjectionParams;
using Source::GridSourceParams;
using Source::RingsSourceParams;
using Source::TessellationSourceParams;
using Source::AuthoredTessellationSourceParams;
using Source::AuthoredProjectedNoiseSourceParams;
using Source::LatticeSourceParams;
using Source::ProjectedNoiseSourceParams;
using Source::SpiralSourceParams;
using Source::TwinWaveSourceParams;
using Source::AuthoredGridSourceParams;
using Source::AuthoredTwinWaveSourceParams;
using Source::AuthoredGridSampleParams;
using Source::AuthoredTwinWaveSampleParams;
using Source::AuthoredRingsSampleParams;
using Source::AuthoredSpiralSampleParams;
using Source::AuthoredLatticeSampleParams;
using Source::AuthoredTessellationSampleParams;
using Source::AuthoredProjectedNoiseSampleParams;
using Surface::DirectSurfaceParams;
using Surface::PeriodicRippleParams;
using Surface::SurfaceNoiseParams;
using Surface::AuthoredSurfaceNoiseParams;
using Surface::AuthoredDirectSurfaceParams;
using Surface::RossbyHaurwitzParams;
using Surface::SymmetricParams;
using Surface::PotentialParams;
using Transfer::IsoValueParams;
using Transfer::SmoothBandsParams;
using Warp::AffineParams;
using Warp::ExplicitAffineParams;
using Warp::PeriodicAffineParams;
using Warp::CurlFlowParams;
using Warp::AuthoredCurlFlowParams;
using Warp::AuthoredWaveShearParams;
using Warp::AuthoredVectorNoiseParams;
using Warp::AuthoredPolarParams;
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
 * @tparam ParamsT Keyed parameter set.
 * @tparam PhaseNoise Whether frames borrow chamber noise patches.
 */
template <typename ParamsT, bool PhaseNoise = false> struct FrameState {
  static constexpr bool PHASE_NOISE_CACHE =
      PhaseNoise; ///< Borrows chamber noise.
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
  ComposedDetail::ResourceStorage<
      ComposedDetail::ResourceFrameBinding<PhaseNoise>::template Type,
      typename ParamsT::ResourceTypes>
      resources; ///< Per-resource clocks and noise views, keyed by resource.
  float palette_oscillation_phase; /**< Phase of the mapping wobble. */
  const PaletteLabLutBase *palette_lab_lut =
      nullptr; ///< Borrowed interval base cache.
  const lenses::DodecahedralTable::Matrix *dodecahedral_matrices =
      nullptr; ///< Runtime-owned immutable reflection matrices.
  const uint8_t *dodecahedral_faces =
      nullptr; ///< Runtime-owned face-major classifiers.
};

/** @brief Column-sampled context with retained color resources.
 * @tparam ParamsT Keyed parameter set.
 * @tparam PhaseNoise Whether frames borrow chamber noise patches. */
template <typename ParamsT, bool PhaseNoise = false>
struct FrameIntervalState : FrameState<ParamsT, PhaseNoise> {
  Color::TemporalPaletteResources temporal; ///< Previous color resources.
};

namespace ComposedDetail {
/** @brief Whether a spec's eligible color modulations sample one phase per
 * column after palette interpolation.
 * @tparam SpecT Effect spec. */
template <typename SpecT>
inline constexpr bool COLOR_PHASE =
    SpecT::TEMPORAL_SAMPLING && SpecT::TEMPORAL_COLOR_PHASE &&
    !SpecT::FUSED_TEMPORAL_COLOR;
/** @brief Whether a spec's palette exposes interval OKLab base caches.
 * @tparam SpecT Effect spec. */
template <typename SpecT>
inline constexpr bool PALETTE_LAB_CACHE =
    SpecT::PALETTE_LAB_CACHE && SpecT::TEMPORAL_SAMPLING &&
    (COLOR_PHASE<SpecT> || SpecT::FUSED_TEMPORAL_COLOR);
} // namespace ComposedDetail

/**
 * @brief Ties a pullback pipeline to one effect's frame state and spec.
 * @details Composed effects render uninstrumented, so the binding pins
 * Pullback::NoInstrumentation. Providers default their spec-selected policies
 * from `Spec`.
 * @tparam FrameT The effect's FrameState specialization.
 * @tparam SpecT The effect's spec.
 */
template <typename FrameT, typename SpecT = Pullback::Spec> struct Binding {
  using FrameState = FrameT;                           ///< Frame type.
  using Instrumentation = Pullback::NoInstrumentation; ///< Uninstrumented.
  using Spec = SpecT;                                  ///< Effect spec.
};

/** @brief Supplies the camera orientation to Pullback::Stage::Rotate.
 * @tparam BindingT Frame binding of the effect. */
template <typename BindingT, typename Family = ProjectionParams>
struct OuterCameraProvider {
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

/** @brief Orientation of one authored sphere-rotation instance.
 * @tparam BindingT Pipeline binding.
 * @tparam Key Rotation resource key. */
template <typename BindingT, ResourceKey Key> struct RotationProvider {
  using Binding = BindingT;                        ///< Pipeline binding.
  using FrameState = typename Binding::FrameState; ///< Bound frame.
  /** @brief Current inverse spin-times-wander orientation.
   * @param frame Current column frame.
   * @return The instance's conjugate orientation. */
  static const math::Quaternion &conjugate(const FrameState &frame) {
    return frame.resources.template get<Key>().conjugate;
  }
};

/**
 * @brief Supplies the projection frame and its parameters to the
 *        Pullback::Projection policies.
 * @tparam BindingT Frame binding of the effect.
 */
template <typename BindingT, typename Family = ProjectionParams>
struct ProjectionProvider {
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
  /** @brief Projection coordinate scale.
   * @param frame Bound frame.
   * @return Authored scale. */
  static float coordinate_scale(const FrameState &frame) {
    return frame.params.template get<"projection">().coordinate_scale;
  }
  /** @brief Peirce strip translation.
   * @param frame Bound frame.
   * @return Authored strip translation. */
  static float layout_scroll(const FrameState &frame) {
    return frame.params.template get<"projection">().layout_scroll;
  }
  /** @brief Bonne standard parallel.
   * @param frame Bound frame.
   * @return Authored parallel, radians. */
  static float standard_parallel(const FrameState &frame) {
    return frame.params.template get<"projection">().standard_parallel;
  }
  /** @brief Polyhedral atlas area correction.
   * @param frame Bound frame.
   * @return Authored blend toward equal area. */
  static float area_correction(const FrameState &frame) {
    return frame.params.template get<"projection">().area_correction;
  }
};

/** @brief Reads an independently keyed lens family.
 * @tparam BindingT Pipeline binding.
 * @tparam Family Lens controls.
 * @tparam Key Parameter resource key. */
template <typename BindingT, typename Family, ResourceKey Key>
struct KeyedLensProvider {
  using Binding = BindingT;                        ///< Pipeline binding.
  using FrameState = typename Binding::FrameState; ///< Bound frame.
  /** @brief Reads authored controls.
   * @param frame Bound frame.
   * @return Keyed controls. */
  static const auto &params(const FrameState &frame) {
    return frame.params.template get<Key>();
  }
};

/** @brief Supplies a Mobius generator and its clock.
 * @tparam BindingT Frame binding.
 * @tparam Key Lens resource key. */
template <typename BindingT, ResourceKey Key = "lens"> struct LensProvider {
  using Binding = BindingT;                        ///< Provider binding.
  using FrameState = typename Binding::FrameState; ///< Bound frame.
  /** @brief Reads authored controls. @param frame Bound frame.
   * @return Keyed generator controls. */
  static const MobiusGeneratorParams &params(const FrameState &frame) {
    return frame.params.template get<Key>();
  }
  /** @brief Reads the animation clock. @param frame Bound frame.
   * @return Keyed clock. */
  static const auto &phase(const FrameState &frame) {
    return frame.resources.template get<Key>();
  }
  /** @brief Reads the interval's shared angle-zero exponential.
   * @param frame Bound frame.
   * @return Shared coefficients, or null for varying controls or magnitude. */
  static const math::MobiusParams *mobius_base(const FrameState &frame) {
    return frame.resources.template get<Key>().mobius_base;
  }
};

/** @brief Supplies the parameters and phase of one animated twist.
 * @tparam BindingT Bound frame.
 * @tparam Key Twist resource key. */
template <typename BindingT, ResourceKey Key> struct TwistProvider {
  using Binding = BindingT;                        ///< Provider binding.
  using FrameState = typename Binding::FrameState; ///< Bound frame.
  /** @brief Reads the twist controls.
   * @param frame Bound frame.
   * @return Keyed parameter family. */
  static const TwistLensParams &params(const FrameState &frame) {
    return frame.params.template get<Key>();
  }
  /** @brief Reads the axis-angle clocks.
   * @param frame Bound frame.
   * @return Keyed phase state. */
  static const Lens::TwistPhase &phase(const FrameState &frame) {
    return frame.resources.template get<Key>();
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
    if constexpr (std::is_same_v<WarpT, PeriodicAffineParams>) {
      const auto windings = Pullback::Warp::cell_windings(params(frame));
      return Pullback::Warp::prepare(
          windings, phase(frame), frame.resources.template get<Key>().rotation,
          windings.lattice_period);
    } else if constexpr (std::is_same_v<WarpT, ExplicitAffineParams>) {
      return Pullback::Warp::prepare(
          params(frame), phase(frame),
          frame.resources.template get<Key>().rotation,
          params(frame).lattice_period);
    } else if constexpr (std::is_same_v<WarpT, AffineParams>) {
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
  using RadialPrepareOwner = SurfaceProvider; ///< Owner of radial preparation.
  using GradientPrepareOwner =
      SurfaceProvider; ///< Owner of gradient qualification.
  /** @brief The surface's noise generator.
   *  @param frame Current frame.
   *  @return The runtime-owned generator. */
  __attribute__((always_inline)) static const FastNoiseLite &
  noise(const FrameState &frame) {
    return *frame.resources.template get<Key>().noise;
  }
  /** @brief Per-frame prepared surface state.
   *  @param frame Current frame.
   *  @return `SphereAnimation::prepare_ripples` for burst ripples,
   *  `Pullback::Surface::prepare_direct` for a directed family,
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
    else {
      auto prepared =
          Pullback::Surface::prepare(frame.resources.template get<Key>().phase);
      if constexpr (std::is_same_v<Family, SurfaceNoiseParams>)
        prepared.unit_default_gradient =
            noise(frame).has_default_gradient_transform();
      return prepared;
    }
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
  using Binding = typename State::Binding;       ///< Frame binding.
  using FrameState = typename State::FrameState; ///< Current frame.
  static constexpr bool LAZY_RADIAL_PREPARE = [] {
    if constexpr (requires {
                    typename State::RadialPrepareOwner;
                    State::radial_tables(std::declval<const FrameState &>());
                    State::radial_fraction(std::declval<const FrameState &>());
                  })
      return std::is_same_v<State, typename State::RadialPrepareOwner>;
    else
      return false;
  }(); ///< Whether table-eligible columns defer exact wave preparation.
  using Prepared = std::conditional_t<
      LAZY_RADIAL_PREPARE, std::optional<SphereAnimation::PreparedRipples>,
      SphereAnimation::PreparedRipples>; ///< Optional exact wave kernels.
  /** @brief Prepares exact waves unless the column has a usable radial table.
   * @param frame Current frame.
   * @return Exact kernels, or an empty optional for deferred preparation. */
  static Prepared prepare(const FrameState &frame) {
    if constexpr (LAZY_RADIAL_PREPARE) {
      const auto *TABLES = State::radial_tables(frame);
      if (TABLES && TABLES->enabled)
        return std::nullopt;
    }
    return State::prepare(frame);
  }
  /** @brief Prepares exact waves for a rejected radial sample.
   * @param input Original input direction.
   * @param frame Current column context.
   * @return Exact grouped displacement and angular travel. */
  HS_HOT_FLASH_MEMBER static SurfaceResult
  exact_fallback(const math::Vector &input, const FrameState &frame) {
    const auto PREPARED = State::prepare(frame);
    return SphereAnimation::grouped_ripples(input, PREPARED,
                                            State::path_length_required(frame));
  }
  /** @brief Applies every active wave in launch order.
   * @param input Unit sphere direction.
   * @param frame Current frame.
   * @param prepared Exact kernels, or an empty optional for radial evaluation.
   * @return Displaced direction and accumulated path length. */
  static SurfaceResult apply(const math::Vector &input, const FrameState &frame,
                             const Prepared &prepared) {
    if constexpr (LAZY_RADIAL_PREPARE) {
      if (prepared)
        return SphereAnimation::grouped_ripples(
            input, *prepared, State::path_length_required(frame));
      return SphereAnimation::radial_ripples(
          input, State::path_length_required(frame),
          State::radial_tables(frame), State::radial_fraction(frame),
          [&](const math::Vector &original)
              HS_HOT_INLINE { return exact_fallback(original, frame); });
    } else if constexpr (
        requires {
          State::radial_tables(frame);
          State::radial_fraction(frame);
        } && !requires { typename State::RadialPrepareOwner; })
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
  __attribute__((always_inline)) static auto prepare(const FrameState &frame) {
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
  /** @brief Number of cosine bands.
   * @param frame Current frame.
   * @return The family's band count. */
  __attribute__((always_inline)) static float
  band_count(const FrameState &frame) {
    return frame.params.template get<Key>().band_count;
  }
  /** @brief Cosine-band phase in radians.
   * @param frame Current frame.
   * @return The family's band phase. */
  __attribute__((always_inline)) static float
  band_phase(const FrameState &frame) {
    return frame.params.template get<Key>().band_phase;
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
 * @tparam FusedV Whether the final temporal color layer may be fused.
 */
template <typename BindingT, Pullback::Color::BrightnessEnvelope BrightnessV,
          bool LabCacheV =
              ComposedDetail::PALETTE_LAB_CACHE<typename BindingT::Spec>,
          bool FusedV = BindingT::Spec::FUSED_TEMPORAL_COLOR>
struct ColorProvider {
  using Binding = BindingT;                        ///< The effect's Binding.
  using FrameState = typename Binding::FrameState; ///< The effect's frame.
  static constexpr bool FUSED_TEMPORAL_COLOR =
      FusedV; ///< Fuses the final layer.
  /** @brief Borrows the interval's palette base cache.
   * @param frame Column frame context.
   * @return Persistent cache, or null for canonical draws. */
  static const PaletteLabLutBase *palette_lab_lut(const FrameState &frame)
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
 * @tparam DynamicV Whether presets select the source and evaluation.
 * @tparam PhaseV Whether eligible sources share column-phase evaluation. */
template <typename BindingT, ResourceKey Key, Pattern::SourceKind SourceV,
          Pattern::Evaluation EvaluationV = Pattern::Evaluation::BAKED_CUBE,
          bool DynamicV = false,
          bool PhaseV = ComposedDetail::COLOR_PHASE<typename BindingT::Spec>>
struct ModulationProvider {
  using Binding = BindingT;                        ///< Pipeline binding.
  using FrameState = typename Binding::FrameState; ///< Frame context.
  static constexpr bool CACHED_NOISE = [] {
    if constexpr (requires { FrameState::PHASE_NOISE_CACHE; })
      return FrameState::PHASE_NOISE_CACHE && !DynamicV &&
             SourceV == Pattern::SourceKind::NOISE &&
             EvaluationV == Pattern::Evaluation::DIRECT;
    return false;
  }(); ///< Whether this source borrows chamber patches.
  /// Whether endpoint noise cubes may share the column color tail.
  static constexpr bool CUBE_PHASE =
      !DynamicV && SourceV == Pattern::SourceKind::NOISE &&
      EvaluationV == Pattern::Evaluation::BAKED_CUBE;
  using PhasePrepared = std::conditional_t<
      CACHED_NOISE, Color::ChamberPhaseModulationState,
      std::conditional_t<CUBE_PHASE, Color::CubePhaseModulationState,
                         Color::PhaseModulationState>>; ///< Column phase state.
  using Prepared = std::conditional_t<PhaseV, PhasePrepared,
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
      PhasePrepared result;
      static_cast<Color::ModulationState &>(result) = BASE;
      if constexpr (requires { frame.temporal.phase_compatible; }) {
        if (!frame.temporal.phase_compatible)
          return result;
        if (!DynamicV &&
            (config.source != SourceV || config.evaluation != EvaluationV))
          return result;
        if (config.mechanic.kind != ColorModulation::Mechanic::HUE_ROTATION &&
            config.mechanic.kind !=
                ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT)
          return result;
        const auto SOURCE = DynamicV ? config.source : SourceV;
        const bool DIRECT = (DynamicV ? config.evaluation : EvaluationV) ==
                            Pattern::Evaluation::DIRECT;
        if (SOURCE == Pattern::SourceKind::DISPLACEMENT)
          result.ready = true;
        else if constexpr (CUBE_PHASE) {
          if (SOURCE == Pattern::SourceKind::NOISE && !DIRECT) {
            if (!resource.modulation || !resource.previous_modulation ||
                !resource.previous_modulation->cube.active ||
                !resource.modulation->cube.active)
              return result;
            result.cube_previous = &resource.previous_modulation->cube;
            result.cube_current = &resource.modulation->cube;
            result.cube_fraction = frame.temporal.fraction;
            result.ready = true;
            return result;
          }
        }
        if (SOURCE == Pattern::SourceKind::NOISE && DIRECT) {
          if (!resource.modulation || !resource.previous_modulation)
            return result;
          if (!Modulation::phase_compatible(
                  resource.previous_modulation->config,
                  resource.modulation->config))
            return result;
          result.noise = &std::get<0>(resource.modulation->source);
          if constexpr (CACHED_NOISE) {
            const auto &a = resource.phase_noise_previous;
            const auto &b = resource.phase_noise_current;
            if (a.patch && b.patch && resource.phase_noise_fraction > 0.0f &&
                resource.phase_noise_fraction < 1.0f &&
                a.patch->matches(*resource.previous_modulation, a.projection) &&
                b.patch->matches(*resource.modulation, b.projection) &&
                Modulation::ChamberNoisePatch::continuous(
                    a.patch->phase, a.projection, b.patch->phase,
                    b.projection)) {
              result.noise_previous = a.patch;
              result.noise_current = b.patch;
              result.noise_fraction = resource.phase_noise_fraction;
              result.noise_to_chamber = frame.projection_conjugate.conjugate();
            }
          }
          const float PHASE =
              Spacetime::ClockInterval(
                  resource.previous_modulation->state.phase.phase_turns,
                  resource.modulation->state.phase.phase_turns, 1.0f)
                  .sample(frame.temporal.fraction);
          if constexpr (CACHED_NOISE)
            if (result.noise_previous && result.noise_current) {
              result.noise_phase = PHASE;
              result.lazy_noise_offset = true;
              result.ready = true;
              return result;
            }
          result.noise_offset = math::noise_sphere_loop_offset(PHASE);
          result.ready = true;
        } else if (SOURCE == Pattern::SourceKind::SPHERE_SMOOTH_GRID &&
                   DIRECT) {
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

/** @brief Color provider with authored palette selection and value envelope.
 * @tparam B Frame binding.
 * @tparam Brightness Fallback envelope for discovery.
 * @tparam LabCache Whether palette interval caches are exposed.
 * @tparam Fused Whether final temporal color may be fused. */
template <typename B,
          Color::BrightnessEnvelope Brightness =
              Color::BrightnessEnvelope::NONE,
          bool LabCache = ComposedDetail::PALETTE_LAB_CACHE<typename B::Spec>,
          bool Fused = B::Spec::FUSED_TEMPORAL_COLOR>
struct AuthoredColorProvider : ColorProvider<B, Brightness, LabCache, Fused> {
  using FrameState = typename B::FrameState;
  /** @brief Current authored brightness envelope.
   * @param frame Current column state.
   * @return Selected envelope. */
  static Color::BrightnessEnvelope
  brightness_envelope(const FrameState &frame) {
    return frame.params.template get<"color">().brightness_envelope;
  }
};
