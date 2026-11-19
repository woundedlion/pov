/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/math/scalar_kernels.h"
#include "core/math/3dmath.h"
#include "core/color/palette_selection.h"
#include "core/render/pullback/contract.h"
#include "core/render/pullback/fields.h"
#include "core/color/hue_rotation_lut.h"
#include "core/color/palette_lab_lut.h"
#include "core/render/pullback/modulation_temporal.h"
#include "core/render/pullback/chamber_noise_cache.h"

/**
 * @file color.h
 * @brief Color and palette policies.
 */

#include "core/render/pullback/color_controls.h"
#include "core/render/pullback/color_state.h"
#include "core/render/pullback/phase_stack.h"
#include "core/render/pullback/palette_cache.h"

namespace Pullback::Color {

/** @brief Generated-palette policy with retained immutable endpoints.
 * @tparam State Palette resource provider. */
template <typename State> struct GeneratedPalette : ApproximationDefaults {
  using SharedBaseOwner = GeneratedPalette<State>;  ///< Cache-contract owner.
  using CachedPhaseOwner = GeneratedPalette<State>; ///< Direct finish owner.
  static constexpr bool PHASE_BLEND_BASE = true;    ///< Common field controls.
  /** @brief Checks the generated-palette provider contract.
   * @tparam Binding Frame binding. */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        {
          State::mapping_weights(frame)
        } -> std::same_as<PaletteMappingWeights>;
        { State::mapping_frequency(frame) } -> std::same_as<float>;
        { State::mapping_phase(frame) } -> std::same_as<float>;
        { State::oscillation_depth(frame) } -> std::same_as<float>;
        { State::oscillation_phase(frame) } -> std::same_as<float>;
        { State::palette(frame) } -> std::same_as<const BakedPalette &>;
        { State::hue_rotation(frame) } -> std::same_as<HueRotationLutView>;
        {
          State::brightness_envelope(frame)
        } -> std::same_as<BrightnessEnvelope>;
        { State::brightness_bottom(frame) } -> std::same_as<float>;
        { State::brightness_top(frame) } -> std::same_as<float>;
        { State::opacity_low(frame) } -> std::same_as<float>;
        { State::opacity_high(frame) } -> std::same_as<float>;
      };
  static constexpr bool TEMPORAL =
      requires(const typename State::Binding::FrameState &frame) {
        State::temporal_resources(frame);
      }; ///< Whether previous palettes are retained.
  static constexpr bool LAB_CACHE =
      requires(const typename State::Binding::FrameState &frame) {
        State::palette_lab_lut(frame);
      }; ///< Whether an RGB-first base cache is available.
  static constexpr bool FUSED_TEMPORAL_COLOR = [] {
    if constexpr (requires { State::FUSED_TEMPORAL_COLOR; })
      return State::FUSED_TEMPORAL_COLOR;
    else
      return false;
  }(); ///< Whether the effect admits one-carrier endpoint color blending.
  static constexpr bool APPROXIMATE =
      LAB_CACHE; ///< Optional base approximation.
  static constexpr ApproximationOracleId ORACLE =
      LAB_CACHE ? ApproximationOracleId::PALETTE_RGB_BLEND_TO_LAB_LUT
                : ApproximationOracleId::NONE; ///< Exact RGB-first reference.
  static constexpr auto METRICS = [] {
    if constexpr (LAB_CACHE)
      return std::array<ApproximationMetric, 3>{
          {{ApproximationDomain::COLOR_CHANNEL,
            ApproximationAggregation::MAXIMUM, 7000.0f, "channel code"},
           {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MEAN,
            256.0f, "channel code"},
           {ApproximationDomain::FRAMEBUFFER, ApproximationAggregation::MAXIMUM,
            4000.0f, "channel code"}}};
    else
      return std::array<ApproximationMetric, 0>{};
  }(); ///< Base-cache qualification bounds.
  using Core = GeneratedPaletteCore<TEMPORAL, LAB_CACHE,
                                    FUSED_TEMPORAL_COLOR>; ///< Sampling core.
  using Prepared = typename Core::Prepared;                ///< Endpoint state.
  /** @brief Binds ordered modulation endpoints against this column's palette.
   * @tparam FrameState Column frame type.
   * @tparam Layers Prepared canonical modulation layer types.
   * @param frame Column frame referencing current and optional retained resources.
   * @param prepared Prepared palette.
   * @param layers Modulation layers in authored order.
   * @return Shared evaluator borrowing the prepared resources. */
  template <LayerLayout... Layout, typename FrameState, typename... Layers>
  static auto bind_modulation_stack(const FrameState &frame,
                                    const Prepared &prepared,
                                    const Layers &...layers)
    requires(sizeof...(Layers) <= MAX_ORDERED_MODULATION_LAYERS &&
             (std::is_base_of_v<ModulationState, Layers> && ...) &&
             (TEMPORAL == requires { frame.temporal.fraction; }) &&
             (!TEMPORAL || requires { frame.temporal.same_palette; }))
  {
    const std::array<const ModulationState *, sizeof...(Layers)> sources{
        &layers...};
    if constexpr (TEMPORAL)
      return ModulationStack<sizeof...(Layers), Layout...>{
          prepared, frame.temporal.same_palette, sources};
    else
      return ModulationStack<sizeof...(Layers), Layout...>{prepared, sources};
  }
  /** @brief Binds modulation layers sampled at the interpolated column phase.
   * @tparam Policies Canonical layer policies in authored order.
   * @param base Prepared palette endpoints.
   * @param fraction Current column weight.
   * @param cached Whether the palette cache supplies the column base.
   * @param layers Ready phase states.
   * @return Scalar-sharing evaluator borrowing the supplied layer states. */
  template <typename... Policies>
  static auto bind_phase_stack(const Prepared &base, float fraction,
                               bool cached,
                               const typename Policies::Prepared &...layers) {
    return PhaseModulationStack{
        base, fraction, cached,
        PhaseLayer<typename Policies::Prepared, Policies::SINE_PHASE>{
            layers}...};
  }
  /** @brief Finishes an eligible cached phase without a deferred carrier.
   * @tparam Layer Canonical phase modulation policy.
   * @param sample Shared field coordinates.
   * @param prepared Prepared palette endpoints.
   * @param modulation Prepared column phase.
   * @param output Receives the final straight-alpha color.
   * @return Whether the RGB cache supplied the final color. */
  template <typename Layer>
  HS_HOT_INLINE static bool
  shade_cached_phase(const FieldSample &sample, const Prepared &prepared,
                     const typename Layer::Prepared &modulation, Color4 &output)
    requires(TEMPORAL && LAB_CACHE &&
             std::is_base_of_v<PhaseModulationState, typename Layer::Prepared>)
  {
    return Color::shade_cached_phase<typename Layer::Prepared,
                                     Layer::SINE_PHASE>(sample, prepared,
                                                        modulation, output);
  }
  /** @brief Prepares immutable endpoint controls.
   * @tparam FrameState Frame context.
   * @param frame Current frame.
   * @return Prepared palette endpoints. */
  template <typename FrameState>
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    GeneratedPaletteState current{
        State::mapping_weights(frame),
        State::mapping_frequency(frame),
        oscillated_mapping_offset(State::mapping_phase(frame),
                                  State::oscillation_depth(frame),
                                  State::oscillation_phase(frame)),
        &State::palette(frame),
        State::hue_rotation(frame),
        State::brightness_envelope(frame),
        State::brightness_bottom(frame),
        State::brightness_top(frame),
        State::opacity_low(frame),
        State::opacity_high(frame)};
    if constexpr (TEMPORAL) {
      const auto resources = State::temporal_resources(frame);
      auto previous = current;
      previous.palette = resources.palette;
      previous.hue_rotation.data = resources.hue_rotation;
      previous.hue_rotation.gamut_mask = resources.hue_rotation_gamut;
      TemporalGeneratedPaletteState result{current, previous,
                                           resources.fraction};
      if constexpr (LAB_CACHE)
        if (const auto *cache = State::palette_lab_lut(frame))
          result.lab = cache->view(resources.fraction);
      return result;
    } else
      return current;
  }
  /** @brief Produces one color endpoint.
   * @tparam FrameState Frame context.
   * @param sample Field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared endpoints.
   * @param previous Selects the retained endpoint.
   * @return Color carrier. */
  template <typename FrameState>
  HS_HOT_INLINE static ColorSample
  apply_endpoint(const FieldSample &sample, const FrameState &frame,
                 const Prepared &prepared, bool previous) {
    (void)frame;
    return Core::apply_endpoint(sample, prepared, previous);
  }
  /** @brief Samples both palettes with their shared mapping and gain controls.
   * @tparam FrameState Frame context.
   * @param sample Shared field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared endpoints.
   * @return Retained endpoint followed by current endpoint. */
  template <typename FrameState>
  HS_HOT_INLINE static std::pair<ColorSample, ColorSample>
  apply_pair(const FieldSample &sample, const FrameState &frame,
             const Prepared &prepared)
    requires TEMPORAL
  {
    (void)frame;
    return Core::apply_pair(sample, prepared);
  }
  /** @brief Tries one cached base without fetching either palette.
   * @tparam FrameState Frame context.
   * @param sample Shared field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared endpoint controls and cache view.
   * @param output Receives the cached color carrier.
   * @return Whether the cache supplied the base. */
  template <typename FrameState>
  HS_HOT_INLINE static bool
  apply_shared_base(const FieldSample &sample, const FrameState &frame,
                    const Prepared &prepared, ColorSample &output)
    requires(TEMPORAL && LAB_CACHE)
  {
    (void)frame;
    return Core::apply_shared_base(sample, prepared, output);
  }
  /** @brief Produces a cached base after column eligibility was checked.
   * @tparam FrameState Frame context.
   * @param sample Shared field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared eligible cache view.
   * @param output Receives the cached color carrier.
   * @pre The base view is active and its fraction is inside (0, 1]. */
  template <typename FrameState>
  HS_HOT_INLINE static void
  apply_shared_base_ready(const FieldSample &sample, const FrameState &frame,
                          const Prepared &prepared, ColorSample &output)
    requires(TEMPORAL && LAB_CACHE)
  {
    (void)frame;
    Core::apply_shared_base_ready(sample, prepared, output);
  }
  /** @brief Blends palette and modulation endpoints into one guarded carrier.
   * @tparam FrameState Frame context.
   * @param sample Shared field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared palette endpoints.
   * @param modulation Prepared modulation endpoints.
   * @param fraction Current endpoint weight.
   * @return Deferred color after canonical endpoint blending. */
  template <bool MixPath = true, typename FrameState>
  HS_HOT_INLINE static ColorSample
  apply_fused_modulation(const FieldSample &sample, const FrameState &frame,
                         const Prepared &prepared,
                         const ModulationState &modulation, float fraction)
    requires(TEMPORAL && FUSED_TEMPORAL_COLOR)
  {
    (void)frame;
    return Core::template apply_fused_modulation<MixPath>(sample, prepared,
                                                          modulation, fraction);
  }
  /** @brief Shades fused endpoints without returning a deferred carrier.
   * @tparam FrameState Frame context.
   * @param sample Shared field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared palette endpoints.
   * @param modulation Prepared modulation endpoints.
   * @param fraction Current endpoint weight.
   * @return Final straight-alpha color. */
  template <bool MixPath = true, typename FrameState>
  HS_HOT_INLINE static Color4
  shade_fused_modulation(const FieldSample &sample, const FrameState &frame,
                         const Prepared &prepared,
                         const ModulationState &modulation, float fraction)
    requires(TEMPORAL && FUSED_TEMPORAL_COLOR)
  {
    (void)frame;
    return Core::template shade_fused_modulation<MixPath>(sample, prepared,
                                                          modulation, fraction);
  }
  /** @brief Produces the current color endpoint.
   * @tparam FrameState Frame context.
   * @param sample Field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared endpoints.
   * @return Color carrier. */
  template <typename FrameState>
  HS_HOT_INLINE static ColorSample apply(const FieldSample &sample,
                                         const FrameState &frame,
                                         const Prepared &prepared) {
    (void)frame;
    return Core::apply(sample, prepared);
  }
};

/** @brief Ordered color modulation policy.
 * @tparam State Layer resource provider. */
template <typename State> struct Modulate : ApproximationDefaults {
  /** @brief Returns the layer layout declared by the provider.
   * @return Static source and mechanic, or the dynamic layout by default. */
  template <typename Owner> static consteval LayerLayout source_layout() {
    if constexpr (requires { State::template source_layout<Owner>(); })
      return State::template source_layout<Owner>();
    return {};
  }
  using PairedOwner = Modulate; ///< Policy admitting paired endpoint updates.
  static constexpr bool APPROXIMATE =
      true; ///< Uses scalar and hue lookup tables.
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS; ///< Reference.
  static constexpr auto METRICS = MODULATION_LUT_METRICS; ///< Error bounds.

  using Prepared = decltype(State::prepare(
      std::declval<
          const typename State::Binding::FrameState &>())); ///< Layer state.
  /// Additional internal approximation qualified separately from the LUT metrics.
  static constexpr bool PHASE_NOISE_APPROXIMATION =
      std::is_same_v<Prepared, ChamberPhaseModulationState>;
  static constexpr bool FUSED_ENDPOINTS =
      std::same_as<Prepared, ModulationState>; ///< Canonical endpoint mode.
  static constexpr bool PATH_MORPH = requires {
    requires State::PATH_MORPH;
  }; ///< Whether a preset edge can blend path and spatial signals.
  static constexpr bool SINE_PHASE = [] {
    if constexpr (requires { State::SINE_PHASE; })
      return State::SINE_PHASE;
    else
      return false;
  }(); ///< Whether the column phase may carry a directional sine.
  /** @brief Checks the modulation provider contract.
   * @tparam Binding Frame binding. */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::prepare(frame) } -> std::same_as<Prepared>;
      };
  /** @brief Prepares source endpoints and continuous controls.
   * @tparam FrameState Frame context.
   * @param frame Frame context.
   * @return Layer endpoints. */
  template <typename FrameState>
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    if constexpr (std::same_as<Prepared, ModulationState>) {
      auto prepared = State::prepare(frame);
      prepared.paired_cube_ready =
          prepared.previous && prepared.current &&
          prepared.previous->evaluate == Modulation::sample_cube &&
          prepared.current->evaluate == Modulation::sample_cube &&
          prepared.previous->cube.active && prepared.current->cube.active;
      return prepared;
    } else
      return State::prepare(frame);
  }
  /** @brief Whether the prepared layer supports shared phase evaluation.
   * @param prepared Layer state.
   * @return Whether this column can use the shared tail. */
  static bool phase_ready(const Prepared &prepared)
    requires std::is_base_of_v<PhaseModulationState, Prepared>
  {
    return prepared.ready;
  }
  /** @brief Applies the column's interpolated analytic source once.
   * @tparam FrameState Frame context.
   * @param sample Carrier to update.
   * @param frame Frame context.
   * @param prepared Layer state.
   * @param scratch Cached analytic samples for this phase. */
  template <typename FrameState>
  HS_HOT_INLINE static void
  apply_phase_in_place(ColorSample &sample, const FrameState &frame,
                       const Prepared &prepared, PhaseScratch &scratch)
    requires std::is_base_of_v<PhaseModulationState, Prepared>
  {
    (void)frame;
    apply_phase_modulation<Prepared, SINE_PHASE>(sample, prepared, scratch);
  }
  /** @brief Applies the column's interpolated analytic source once.
   * @tparam FrameState Frame context.
   * @param sample Carrier to update.
   * @param frame Frame context.
   * @param prepared Layer state. */
  template <typename FrameState>
  HS_HOT_INLINE static void apply_phase_in_place(ColorSample &sample,
                                                 const FrameState &frame,
                                                 const Prepared &prepared)
    requires std::is_base_of_v<PhaseModulationState, Prepared>
  {
    PhaseScratch scratch;
    apply_phase_in_place(sample, frame, prepared, scratch);
  }
  /** @brief Samples a prepared column's analytic phase source.
   * @param field Source coordinates and displacement travel.
   * @param prepared Ready column phase state.
   * @return Source scalar before the mechanic is applied. */
  HS_HOT_INLINE static float sample_phase_scalar(const FieldSample &field,
                                                 const Prepared &prepared)
    requires std::is_base_of_v<PhaseModulationState, Prepared>
  {
    PhaseScratch scratch;
    return sample_phase_scalar(field, prepared, scratch);
  }
  /** @brief Samples a prepared column's analytic phase source.
   * @param field Source coordinates and displacement travel.
   * @param prepared Ready column phase state.
   * @param scratch Lookups shared with the tail's other layers.
   * @return Source scalar before the mechanic is applied. */
  HS_HOT_INLINE static float sample_phase_scalar(const FieldSample &field,
                                                 const Prepared &prepared,
                                                 PhaseScratch &scratch)
    requires std::is_base_of_v<PhaseModulationState, Prepared>
  {
    return Color::sample_phase_scalar<Prepared, SINE_PHASE>(field, prepared,
                                                            scratch);
  }
  /** @brief Modulates one endpoint.
   * @tparam FrameState Frame context.
   * @param input Color carrier.
   * @param frame Frame context.
   * @param prepared Layer endpoints.
   * @param previous Selects the retained endpoint.
   * @return Modulated carrier. */
  template <typename FrameState>
  static ColorSample apply_endpoint(const ColorSample &input,
                                    const FrameState &frame,
                                    const Prepared &prepared, bool previous) {
    (void)frame;
    return apply_modulation<PATH_MORPH>(input, prepared, previous);
  }
  /** @brief Modulates one endpoint without copying its carrier.
   * @tparam FrameState Frame context.
   * @param sample Color carrier to update.
   * @param frame Frame context.
   * @param prepared Layer endpoints.
   * @param previous Selects the retained endpoint. */
  template <typename FrameState>
  static void apply_endpoint_in_place(ColorSample &sample,
                                      const FrameState &frame,
                                      const Prepared &prepared, bool previous) {
    (void)frame;
    apply_modulation_in_place<PATH_MORPH>(sample, prepared, previous);
  }

  /** @brief Updates both endpoints through one prepared layer.
   * @tparam FrameState Frame context.
   * @param previous Retained carrier to update.
   * @param current Current carrier to update.
   * @param frame Frame context.
   * @param prepared Prepared layer state. */
  template <typename FrameState>
  static void
  apply_endpoints_in_place(ColorSample &previous, ColorSample &current,
                           const FrameState &frame, const Prepared &prepared) {
    (void)frame;
    apply_modulation_pair_in_place<PATH_MORPH>(previous, current, prepared);
  }
  /** @brief Modulates the current endpoint.
   * @tparam FrameState Frame context.
   * @param input Color carrier.
   * @param frame Frame context.
   * @param prepared Layer endpoints.
   * @return Modulated carrier. */
  template <typename FrameState>
  static ColorSample apply(const ColorSample &input, const FrameState &frame,
                           const Prepared &prepared) {
    return apply_endpoint(input, frame, prepared, false);
  }
};

} // namespace Pullback::Color
