/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstddef>
#include <type_traits>
#include "core/color/authored_palette.h"
#include "core/color/gamut_lut.h"
#include "core/color/palette_cycler.h"
#include "core/control/params.h"
#include "core/memory.h"
#include "core/render/pullback/composed_descriptors.h"
#include "core/render/pullback/composed_interval.h"
#include "core/render/pullback/composed_spec_policies.h"

/** @file composed_traits.h
 * @brief `ComposedTraits`, the compile-time policy of a composed effect.
 */

namespace Pullback {

/**
 * @brief Policy constants, derived types and spec rules of one composed
 *        effect; holds no state.
 * @tparam W Canvas width in pixels.
 * @tparam SpecT The effect's `Spec`.
 */
template <int W, typename SpecT> struct ComposedTraits {
  using Params = ParamsFor<SpecT>; ///< Keyed parameter storage of the Spec.
  using ColorFamily = std::remove_cvref_t<
      decltype(std::declval<Params &>()
                   .template get<"color">())>; ///< Color controls.
  static constexpr bool AUTHORED_PALETTE =
      std::is_same_v<ColorFamily,
                     Color::AuthoredColorParams>; ///< Full authored sources.
  using Spec = SpecT;                  ///< The effect's ranked stage Spec.
  using Metadata = MetadataFor<SpecT>; ///< Topology derived from typed stages.
  /// Whether runtime identities follow authored nodes.
  static constexpr bool AUTHORED_TOPOLOGY = SpecT::AUTHORED_TOPOLOGY;
  /// Whether the effect samples time across each two-arm half-turn.
  static constexpr bool TEMPORAL_SAMPLING = SpecT::TEMPORAL.sampling();
  /// Whether eligible color modulations sample one phase per column.
  static constexpr bool COLOR_PHASE = ComposedDetail::COLOR_PHASE<SpecT>;
  /// Whether an eligible final temporal color layer is fused.
  static constexpr bool FUSED_TEMPORAL_COLOR =
      ComposedDetail::FUSED_TEMPORAL_COLOR<SpecT>;
  /// Whether the palette retains interval OKLab bases.
  static constexpr bool PALETTE_LAB_CACHE =
      ComposedDetail::PALETTE_LAB_CACHE<SpecT>;
  static_assert(Metadata::ColorModes::LAB_CACHE == PALETTE_LAB_CACHE,
                "palette base cache policy must match its runtime storage");
  static_assert(Metadata::ColorModes::FUSED_TEMPORAL_COLOR ==
                    FUSED_TEMPORAL_COLOR,
                "temporal color fusion policy must match its spec");
  /// Whether the palette retains final hue colors over its interval bases.
  static constexpr bool PALETTE_RGB_CACHE =
      SpecT::TEMPORAL.cache() == PaletteCache::LAB_RGB ||
      SpecT::TEMPORAL.cache() == PaletteCache::LAB_RGB_MIDPOINT;
  /// Whether the final hue cache stores only one admitted midpoint row.
  static constexpr bool PALETTE_RGB_MIDPOINT_ONLY =
      SpecT::TEMPORAL.cache() == PaletteCache::LAB_RGB_MIDPOINT;
  static_assert(!PALETTE_RGB_CACHE ||
                    (SpecT::PALETTE_LAB_VALUES ==
                         HueRotationLutView::VALUE_STEPS &&
                     Metadata::BRIGHTNESS == Color::BrightnessEnvelope::NONE),
                "final hue caches require the hue-rotation value grid and no "
                "brightness envelope");
  /// Whether chamber direct noise is cached at temporal endpoints.
  static constexpr bool PHASE_NOISE_CACHE =
      SpecT::TEMPORAL.noise() == PhaseNoise::CHAMBER_CACHE;
  static_assert(!PHASE_NOISE_CACHE ||
                    std::is_base_of_v<Lens::DodecahedralKaleidoscope,
                                      typename Metadata::LensPolicy>,
                "phase noise cache requires a dodecahedral chamber");
  /// Whether the dodecahedral reflection matrices are copied into RAM.
  static constexpr bool CACHE_DODECAHEDRAL_MATRICES =
      SpecT::DODECAHEDRAL_CACHE.matrices;
  /// Whether the dodecahedral cube classifiers are copied into RAM.
  static constexpr bool CACHE_DODECAHEDRAL_FACES =
      SpecT::DODECAHEDRAL_CACHE.faces;
  static_assert(!CACHE_DODECAHEDRAL_MATRICES ||
                    std::is_base_of_v<Lens::DodecahedralTableKaleidoscope,
                                      typename Metadata::LensPolicy>,
                "cached dodecahedral matrices require a dodecahedral lens");
  static_assert(!CACHE_DODECAHEDRAL_FACES ||
                    std::is_base_of_v<Lens::DodecahedralTableKaleidoscope,
                                      typename Metadata::LensPolicy>,
                "cached cube classifiers require a dodecahedral lens");
  static_assert(!TEMPORAL_SAMPLING || (W >= 2 && W % 2 == 0),
                "Column-time sampling requires an even two-arm canvas width");
  using SourceFamily = typename Params::template Family<"source">;
  /// Whether the source policy blends authored endpoint fields.
  static constexpr bool SOURCE_MORPH =
      SpecT::SOURCE_MORPH || requires { typename SourceFamily::PlanarWeights; };
  using SpatialFrame =
      std::conditional_t<TEMPORAL_SAMPLING,
                         FrameIntervalState<Params, PHASE_NOISE_CACHE>,
                         Pullback::FrameState<Params, PHASE_NOISE_CACHE>>;
  /// Frame context, including retained color resources when sampling time.
  using FrameState = std::conditional_t<
      SpecT::HEX_MORPH_TABLE, Source::HexMorphFrame<SpatialFrame>,
      std::conditional_t<SOURCE_MORPH, Source::EndpointMorphFrame<SpatialFrame>,
                         SpatialFrame>>;

  static_assert(
      !SpecT::HEX_MORPH_TABLE ||
      (SOURCE_MORPH &&
       std::is_base_of_v<
           Source::FixedTessellationParams<Pattern::TessellationKind::HEXAGONAL,
                                           Pattern::TessellationFill::CELLS>,
           SourceFamily>));
  using Binding = Pullback::Binding<FrameState, SpecT>; ///< Pipeline binding.
  /// Whether projection controls are present.
  static constexpr bool HAS_FRAME_PARAMS = Params::template HAS<"projection">;
  /// Whether the projection orientation follows a random walk and spin.
  static constexpr bool ANIMATED_PROJECTION =
      SpecT::PROJECTION_MOTION != ProjectionMotion::STATIC && HAS_FRAME_PARAMS;
  /// Whether the projection frame spins and wanders on the tilted base.
  static constexpr bool TILTED_PROJECTION_FRAME =
      SpecT::PROJECTION_MOTION == ProjectionMotion::ANIMATED_TILTED;
  static_assert(!TILTED_PROJECTION_FRAME ||
                    (ANIMATED_PROJECTION &&
                     !requires(Params p) {
                       p.template get<"projection">().frame;
                     }),
                "ANIMATED_TILTED motion needs projection controls without an "
                "authored frame parameter");
  /// Whether the effect retains palette hue tables.
  static constexpr bool HUE_ROTATION_LUT =
      Metadata::MODULATION && SpecT::HUE_ROTATION_LUT;
  /**
   * @brief Whether any resource of a kind samples a noise field.
   * @tparam Kind Resource kind to test.
   * @return True when a `Kind` resource satisfies
   *         `ComposedDetail::RESOURCE_NOISE`.
   */
  template <ResourceKind Kind> static consteval bool has_noise() {
    bool result = false;
    Params{}.visit([&]<typename Resource>(const auto &) {
      if constexpr (Resource::KIND == Kind &&
                    ComposedDetail::RESOURCE_NOISE<Resource>)
        result = true;
    });
    return result;
  }
  /// Whether a surface resource samples noise.
  static constexpr bool HAS_SURFACE_NOISE = has_noise<ResourceKind::SURFACE>();
  /// Whether a warp resource samples noise.
  static constexpr bool HAS_WARP_NOISE = has_noise<ResourceKind::WARP>();
  /// Whether a source resource samples noise.
  static constexpr bool HAS_SOURCE_NOISE = has_noise<ResourceKind::SOURCE>();
  /// The Spec's pipeline bound to this effect's frame.
  using RenderPipeline = typename SpecT::template Pipeline<Binding>;
  static_assert(!SpecT::FACTORED_COLUMNS || RenderPipeline::LINEAR_COLUMNS,
                "factored columns require a rotation/projection prefix");
  static_assert((Metadata::COVERAGE != ProjectionCoverageMode::EDGE_FADE &&
                 Metadata::COVERAGE !=
                     ProjectionCoverageMode::EDGE_FADE_SQUARED) ||
                    Metadata::EDGE_DISTANCE_AVAILABLE,
                "edge-fade coverage requires projection edge distance");
  static_assert(!Metadata::DISPLACEMENT || Metadata::PATH_TRACKED,
                "displacement modulation requires pipeline path tracking");
  /// Prepared per-frame state the pipeline shades from.
  using Frame = typename RenderPipeline::Frame;

  /** Descriptors the arena-backed parameter array holds; every slider an effect
      registers has to fit. */
  static constexpr size_t PARAM_CAPACITY =
      ComposedDetail::parameter_count<SpecT>(typename Params::ResourceTypes{});

  /**
   * @brief Persistent-arena bytes one effect instance claims.
   * @details init() takes the palette cycler's generated arena, the external
   * ParamDef array, one State, qualified parameter names and, on device, the
   * full gamut LUT. An authored-palette effect's first draw bakes one palette;
   * a temporal effect's clones one palette bake.
   * @tparam State The runtime's persistent aggregate.
   * @return Upper bound on the instance's persistent allocations.
   */
  template <typename State> static consteval size_t footprint_bytes() {
    return (SpecT::CACHE_FULL_GAMUT_LUT
                ? gamut_lut_bytes(GAMUT_LUT_ANGLE_STEPS, GAMUT_LUT_L_STEPS)
                : 0) +
           PaletteCycler::generated_arena_bytes() +
           (AUTHORED_PALETTE ? BakedPalette::required_arena_bytes() : 0) +
           PARAM_CAPACITY * sizeof(ParamDef) + sizeof(State) + alignof(State) +
           (TEMPORAL_SAMPLING ? BakedPalette::required_arena_bytes() : 0) +
           alignof(ParamDef) +
           ComposedDetail::parameter_name_bytes<SpecT>(
               typename Params::ResourceTypes{});
  }
};

} // namespace Pullback
