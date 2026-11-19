/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Prepared palette and modulation endpoint state. */

#include "core/render/pullback/color_controls.h"
#include "core/math/scalar_kernels.h"
#include "core/math/3dmath.h"
#include "core/color/palette_selection.h"
#include "core/render/pullback/contract.h"
#include "core/render/pullback/fields.h"
#include "core/color/hue_rotation_lut.h"
#include "core/color/palette_lab_lut.h"
#include "core/render/pullback/modulation_temporal.h"
#include "core/render/pullback/chamber_noise_cache.h"

namespace Pullback::Color {

/** @brief Per-frame generated-palette controls and immutable resources. */
struct GeneratedPaletteState {
  PaletteMappingWeights mapping;   ///< Mapping weights.
  float mapping_frequency;         ///< Palette repeats.
  float mapping_offset;            ///< Mapping phase.
  const BakedPalette *palette;     ///< Base palette.
  HueRotationLutView hue_rotation; ///< Optional base-palette hue table.
  BrightnessEnvelope brightness_envelope =
      BrightnessEnvelope::NONE;   ///< Gain shape.
  float brightness_bottom = 0.0f; ///< Low gain.
  float brightness_top = 1.0f;    ///< High gain.
  float opacity_low = 1.0f;       ///< Alpha at value zero.
  float opacity_high = 1.0f;      ///< Alpha at value one.
};
/** @brief Retained palette resources and column progress. */
struct TemporalPaletteResources {
  const BakedPalette *palette = nullptr; ///< Previous palette.
  const Pixel *hue_rotation = nullptr;   ///< Previous hue table.
  float fraction = 1.0f;                 ///< Current endpoint weight.
  bool same_palette = false; ///< Identical endpoint palette entries.
  bool phase_compatible =
      false; ///< Stable modulation recipes across the interval.
  const HueRotationGamutMask *hue_rotation_gamut =
      nullptr; ///< Previous hue cell flags.
};
/** @brief Prepared generated-palette endpoint pair. */
struct TemporalGeneratedPaletteState {
  GeneratedPaletteState current;  ///< Current palette.
  GeneratedPaletteState previous; ///< Retained palette.
  float fraction = 1.0f;          ///< Current endpoint weight.
  PaletteLabTimeView lab{};       ///< Optional interval base cache.
};
/** @brief Samples an endpoint with shared palette controls.
 * @param sample Shared field coordinates.
 * @param state Endpoint resources.
 * @param value Prepared palette coordinate.
 * @param opacity Shared coverage and opacity.
 * @param gain Shared final brightness.
 * @return Canonical deferred endpoint. */
HS_HOT_INLINE inline ColorSample
apply_palette_endpoint_impl(const FieldSample &sample,
                            const GeneratedPaletteState &state, float value,
                            float opacity, float gain) {
  Color4 color = state.palette->get(value);
  color.alpha *= opacity;
  ColorSample result{
      color,
      sample,
      state.hue_rotation.active ? state.hue_rotation.data : nullptr,
      value,
      0.0f,
      gain};
  result.hue_rotation_gamut = state.hue_rotation.gamut_mask;
  return result;
}
/** @brief Samples a canonical endpoint from ITCM.
 * @param sample Shared field coordinates.
 * @param state Endpoint resources.
 * @param value Prepared palette coordinate.
 * @param opacity Shared coverage and opacity.
 * @param gain Shared final brightness.
 * @return Canonical deferred endpoint. */
HS_HOT_ITCM_MEMBER inline ColorSample
apply_palette_endpoint(const FieldSample &sample,
                       const GeneratedPaletteState &state, float value,
                       float opacity, float gain) {
  return apply_palette_endpoint_impl(sample, state, value, opacity, gain);
}
/** @brief Resolves palette opacity at the sampled field value.
 * @param sample Sampled field and coverage.
 * @param state Palette opacity controls.
 * @return Coverage-scaled opacity. */
HS_HOT_INLINE inline float palette_opacity(const FieldSample &sample,
                                           const GeneratedPaletteState &state) {
  return sample.coverage *
         hs::lerp(state.opacity_low, state.opacity_high, sample.value);
}
/** @brief Maps a sampled field value onto the palette coordinate.
 * @param sample Sampled field.
 * @param state Palette mapping controls.
 * @return Prepared palette coordinate. */
HS_HOT_INLINE inline float palette_value(const FieldSample &sample,
                                         const GeneratedPaletteState &state) {
  return palette_mapping_coordinate(sample.value, state.mapping,
                                    state.mapping_frequency,
                                    state.mapping_offset);
}
/** @brief Resolves the envelope brightness gain at the sampled field value.
 * @param sample Sampled field.
 * @param state Palette envelope controls.
 * @return Final brightness gain. */
HS_HOT_INLINE inline float palette_gain(const FieldSample &sample,
                                        const GeneratedPaletteState &state) {
  return brightness_envelope_gain(sample.value, state.brightness_envelope,
                                  state.brightness_bottom,
                                  state.brightness_top);
}
/** @brief Palette controls shared by every endpoint of one sample. */
struct PaletteSampleControls {
  float value;   ///< Prepared palette coordinate.
  float opacity; ///< Coverage-scaled opacity.
  float gain;    ///< Final brightness gain.
};
/** @brief Resolves the palette coordinate, opacity and gain of a sample.
 * @param sample Sampled field and coverage.
 * @param state Palette controls.
 * @return Shared endpoint controls. */
HS_HOT_INLINE inline PaletteSampleControls
palette_controls(const FieldSample &sample,
                 const GeneratedPaletteState &state) {
  return {palette_value(sample, state), palette_opacity(sample, state),
          palette_gain(sample, state)};
}
/** @brief Maps a sampled field into the color carrier.
 * @param sample Sampled field.
 * @param state Palette controls and resources.
 * @return Base color with deferred brightness. */
HS_HOT_ITCM_MEMBER inline ColorSample
apply_generated_palette(const FieldSample &sample,
                        const GeneratedPaletteState &state) {
  const auto controls = palette_controls(sample, state);
  return apply_palette_endpoint_impl(sample, state, controls.value,
                                     controls.opacity, controls.gain);
}
/** @brief Completes deferred color adjustments.
 * @param sample Color carrier.
 * @return Final straight-alpha color. */
inline Color4 finish(const ColorSample &sample) {
  return finish_color_sample(sample);
}

/** @brief Blends two authored mechanic contributions into one unit-strength scalar.
 * @param from_scalar Outgoing endpoint's signal.
 * @param from_amount Outgoing mechanic strength.
 * @param to_scalar Incoming endpoint's signal.
 * @param to_amount Incoming mechanic strength.
 * @param fraction Incoming endpoint weight.
 * @return Equivalent signal for a unit-strength mechanic. */
template <ColorModulation::Mechanic Kind>
HS_HOT_INLINE inline float
blend_modulation_scalar(float from_scalar, float from_amount, float to_scalar,
                        float to_amount, float fraction) {
  using M = ColorModulation::Mechanic;
  if constexpr (Kind == M::HUE_ROTATION)
    return math::wrap_t(
        hs::lerp(from_amount * from_scalar, to_amount * to_scalar, fraction));
  else if constexpr (Kind == M::POSITIVE_LIGHTNESS_LIFT)
    return hs::lerp(from_amount * hs::clamp(from_scalar, 0.0f, 1.0f),
                    to_amount * hs::clamp(to_scalar, 0.0f, 1.0f), fraction);
  else {
    static_assert(Kind == M::OPACITY);
    const auto opacity = [](float scalar, float amount) {
      return ColorModulation::Opacity::apply(Color4(Pixel(0, 0, 0), 1.0f),
                                             scalar, {amount})
          .alpha;
    };
    return 2.0f * hs::lerp(opacity(from_scalar, from_amount),
                           opacity(to_scalar, to_amount), fraction) -
           1.0f;
  }
}

/** @brief Blends contributions for a runtime mechanic.
 * @param kind Color transformation.
 * @param from_scalar Outgoing endpoint's signal.
 * @param from_amount Outgoing mechanic strength.
 * @param to_scalar Incoming endpoint's signal.
 * @param to_amount Incoming mechanic strength.
 * @param fraction Incoming endpoint weight.
 * @return Equivalent signal for a unit-strength mechanic. */
HS_HOT_FLASH_MEMBER inline float
blend_modulation_scalar(ColorModulation::Mechanic kind, float from_scalar,
                        float from_amount, float to_scalar, float to_amount,
                        float fraction) {
  using M = ColorModulation::Mechanic;
  if (kind == M::HUE_ROTATION)
    return blend_modulation_scalar<M::HUE_ROTATION>(
        from_scalar, from_amount, to_scalar, to_amount, fraction);
  if (kind == M::POSITIVE_LIGHTNESS_LIFT)
    return blend_modulation_scalar<M::POSITIVE_LIGHTNESS_LIFT>(
        from_scalar, from_amount, to_scalar, to_amount, fraction);
  return blend_modulation_scalar<M::OPACITY>(from_scalar, from_amount,
                                             to_scalar, to_amount, fraction);
}

/** @brief Authored spatial and displacement signals for one changing layer. */
struct PathModulationMix {
  /// Mechanic strength of the spatial endpoint.
  float spatial_amount = 0.0f;
  /// Mechanic strength of the displacement endpoint.
  float path_amount = 0.0f;
  /// Incoming endpoint weight, or negative when inactive.
  float fraction = -1.0f;
  /// Whether the displacement endpoint is outgoing.
  bool path_first = false;

  /// @return Whether the layer carries a source transition.
  bool active() const { return fraction >= 0.0f; }

  /** @brief Blends source signals for a statically selected mechanic.
   * @param spatial Spatial source signal.
   * @param path Displacement path length.
   * @return Equivalent signal for a unit-strength mechanic. */
  template <ColorModulation::Mechanic Kind>
  HS_HOT_INLINE float sample(float spatial, float path) const {
    return sample_with(blend_modulation_scalar<Kind>, spatial, path);
  }

  /** @brief Blends source signals for a runtime mechanic.
   * @param kind Color transformation.
   * @param spatial Spatial source signal.
   * @param path Displacement path length.
   * @return Equivalent signal for a unit-strength mechanic. */
  HS_HOT_INLINE float sample(ColorModulation::Mechanic kind, float spatial,
                             float path) const {
    return sample_with(
        [kind](float a, float amount_a, float b, float amount_b, float t) {
          return blend_modulation_scalar(kind, a, amount_a, b, amount_b, t);
        },
        spatial, path);
  }

private:
  template <typename Blend>
  HS_HOT_INLINE float sample_with(Blend blend, float spatial,
                                  float path) const {
    if (path_first)
      return blend(path, path_amount, spatial, spatial_amount, fraction);
    return blend(spatial, spatial_amount, path, path_amount, fraction);
  }
};

/** @brief One color modulation's prepared endpoint pair. */
struct ModulationState {
  const Modulation::Prepared *current = nullptr;  ///< Current spatial source.
  const Modulation::Prepared *previous = nullptr; ///< Retained spatial source.
  Pattern::SourceKind source = Pattern::SourceKind::NOISE; ///< Scalar source.
  ColorModulation::Params mechanic; ///< Color operation and strength.
  bool enabled = false;             ///< Whether this layer contributes.
  bool paired_cube_ready =
      false; ///< Immutable endpoints admit shared cube taps.
  /// Transition between a spatial and a displacement signal.
  PathModulationMix path_mix{};
  /** @brief Selects or blends the signal this layer consumes.
   * @param spatial Spatial source signal.
   * @param path Displacement path length.
   * @return Signal after any enabled source transition. */
  template <bool MixPath = true>
  HS_HOT_INLINE float scalar(float spatial, float path) const {
    if constexpr (MixPath)
      if (path_mix.active())
        return path_mix.sample(mechanic.kind, spatial, path);
    return spatial;
  }
};

/** @brief Interval-base cache tier of a borrowed palette view. */
enum class CacheTier : uint8_t {
  NONE, ///< No OKLab rows, or only one RGB hue table.
  LAB,  ///< OKLab rows without RGB hue tables.
  RGB   ///< OKLab rows with both RGB hue tables.
};

/** @brief Cache tier of a borrowed palette view.
 * @param lab Interval cache view.
 * @return The view's tier. */
HS_HOT_INLINE inline CacheTier cache_tier(const PaletteLabTimeView &lab) {
  if (!lab.low)
    return CacheTier::NONE;
  if (lab.rgb_low && lab.rgb_high)
    return CacheTier::RGB;
  return !lab.rgb_low && !lab.rgb_high ? CacheTier::LAB : CacheTier::NONE;
}

/** @brief Whether a modulation layer changes color.
 * @param modulation Layer controls.
 * @return Whether the layer is enabled with nonzero strength. */
HS_HOT_INLINE inline bool modulation_active(const ModulationState &modulation) {
  return ColorModulation::active(modulation.enabled, modulation.mechanic);
}

/** @brief Whether a modulation layer is inert or only rotates hue.
 * @param modulation Layer controls.
 * @return Whether the layer leaves everything but hue unchanged. */
HS_HOT_INLINE inline bool hue_only(const ModulationState &modulation) {
  return ColorModulation::hue_only(modulation.enabled, modulation.mechanic);
}

/** @brief Finishing route of a cached palette followed by one layer. */
enum class CachedColorPath : uint8_t {
  CANONICAL,       ///< Canonical endpoint blending.
  RGB_PHASE,       ///< RGB hue tables at one column phase.
  LAB_ENDPOINT,    ///< OKLab rows at the current endpoint.
  LAB_HUE_INTERVAL ///< OKLab rows with an interpolated hue.
};

/** @brief Selects the cached finishing route of a palette and one layer.
 * @param prepared Palette endpoints and cache view.
 * @param modulation Layer controls.
 * @param fraction Column endpoint weight.
 * @param phase_ready Whether the layer admits one shared column phase.
 * @return The admitted route, or CANONICAL. */
HS_HOT_INLINE inline CachedColorPath
cached_color_path(const TemporalGeneratedPaletteState &prepared,
                  const ModulationState &modulation, float fraction,
                  bool phase_ready) {
  const bool flat =
      prepared.current.brightness_envelope == BrightnessEnvelope::NONE;
  switch (cache_tier(prepared.lab)) {
  case CacheTier::RGB:
    return interval_fraction(prepared.fraction) && flat && phase_ready &&
                   hue_only(modulation)
               ? CachedColorPath::RGB_PHASE
               : CachedColorPath::CANONICAL;
  case CacheTier::LAB:
    if (fraction == 1.0f && prepared.fraction == 1.0f && flat &&
        hue_only(modulation))
      return CachedColorPath::LAB_ENDPOINT;
    if (fraction == prepared.fraction && fraction > 0.0f && fraction < 1.0f &&
        modulation_active(modulation) &&
        modulation.mechanic.kind == ColorModulation::Mechanic::HUE_ROTATION &&
        flat)
      return CachedColorPath::LAB_HUE_INTERVAL;
    return CachedColorPath::CANONICAL;
  case CacheTier::NONE:
    break;
  }
  return CachedColorPath::CANONICAL;
}

/** @brief Eligible spatial source prepared at the column phase. */
struct PhaseModulationState : ModulationState {
  /// Interpolated axial phases; a directional sine uses the first wave.
  Pattern::SphereSmoothGrid::Prepared grid{};
  const Pattern::SphereNoise::Prepared *noise =
      nullptr; ///< Direct generator; stationary tables borrow current->immutable_table.
  math::Vector noise_offset{}; ///< Noise loop at the column phase.
  bool ready = false;          ///< Whether one shared color tail is admitted.
};

/** @brief Endpoint noise cubes blended by a baked-noise color tail. */
struct CubePhaseModulationState : PhaseModulationState {
  const Pattern::SphereLutView *cube_previous = nullptr; ///< Retained cube.
  const Pattern::SphereLutView *cube_current = nullptr;  ///< Current cube.
  float cube_fraction = 0.0f; ///< Current cube weight.
};
static_assert(sizeof(Pattern::SphereSmoothGrid::Prepared) == 40);

/** @brief Optional chamber patches shared by a direct-noise color tail. */
struct ChamberPhaseModulationState : PhaseModulationState {
  const Modulation::ChamberNoisePatch *noise_previous =
      nullptr; ///< Retained samples.
  const Modulation::ChamberNoisePatch *noise_current =
      nullptr;                       ///< Current samples.
  math::Quaternion noise_to_chamber; ///< Inverse column projection.
  float noise_fraction = 0.0f;       ///< Current endpoint weight.
  float noise_phase = 0.0f;          ///< Loop phase for uncached lookups.
  bool lazy_noise_offset = false;    ///< Resolve the loop offset on fallback.
};

/** @brief Samples one prepared phase, sharing per-pixel cube coordinates.
 * @tparam Phase Prepared modulation phase type.
 * @tparam SinePhase Whether to sample directional-sine phase sources.
 * @param field Source coordinates and displacement travel.
 * @param prepared Ready column phase.
 * @param scratch Lookups shared by this phase tail.
 * @return Source scalar before its mechanic. */
template <typename Phase, bool SinePhase = true>
HS_HOT_INLINE inline float sample_phase_scalar(const FieldSample &field,
                                               const Phase &prepared,
                                               PhaseScratch &scratch)
  requires std::is_base_of_v<PhaseModulationState, Phase>
{
  float scalar = 0.0f;
  if (prepared.source == Pattern::SourceKind::DISPLACEMENT)
    scalar = field.path_length;
  else if (prepared.source == Pattern::SourceKind::NOISE) {
    if constexpr (std::is_same_v<Phase, CubePhaseModulationState>)
      if (prepared.cube_current) {
        const auto &cube = scratch.tap(field.sphere);
        return hs::lerp(
            Pattern::bounded_scalar(cube.sample(*prepared.cube_previous)),
            Pattern::bounded_scalar(cube.sample(*prepared.cube_current)),
            prepared.cube_fraction);
      }
    if (!prepared.noise)
      return Pattern::bounded_scalar(
          prepared.current->immutable_table->sample_tap(
              *prepared.current, scratch.tap(field.sphere)));
    bool cached = false;
    if constexpr (std::is_same_v<Phase, ChamberPhaseModulationState>)
      if (prepared.noise_previous && prepared.noise_current)
        cached = Modulation::sample_chamber_noise(
            field.sphere, prepared.noise_to_chamber, *prepared.noise_previous,
            *prepared.noise_current, prepared.noise_fraction, scalar);
    if (!cached) {
      const auto &noise = *prepared.noise;
      auto noise_offset = prepared.noise_offset;
      if constexpr (std::is_same_v<Phase, ChamberPhaseModulationState>)
        if (prepared.lazy_noise_offset)
          noise_offset = math::noise_sphere_loop_offset(prepared.noise_phase);
      const auto q = math::noise_sphere_coordinate(
          scratch.admit(field.sphere).value(), noise.scale, noise_offset);
      scalar =
          Pattern::bounded_scalar(noise.noise.GetNoiseSingle(q.x, q.y, q.z));
    }
  } else {
    if constexpr (SinePhase)
      if (prepared.source == Pattern::SourceKind::DIRECTIONAL_SINE)
        return Pattern::DirectionalSine::sample(scratch.admit(field.sphere),
                                                prepared.grid.a);
    scalar = Pattern::SphereSmoothGrid::sample(scratch.admit(field.sphere),
                                               prepared.grid);
  }
  return scalar;
}

/** @brief Accumulates one deferred hue rotation.
 * @param hue Prior hue in turns.
 * @param scalar Source value.
 * @param mechanic Rotation controls.
 * @return Wrapped hue in turns. */
HS_HOT_INLINE inline float
accumulate_hue(float hue, float scalar,
               const ColorModulation::Params &mechanic) {
  return math::wrap_t(hue + math::wrap_t(mechanic.amount * scalar));
}

/** @brief Samples a cached RGB palette with the accumulated phase hue.
 * @param sample Field coordinates and coverage.
 * @param prepared Cached interval base and palette controls.
 * @param hue Accumulated hue in turns.
 * @return Final straight-alpha color. */
HS_HOT_INLINE inline Color4
cached_phase_color(const FieldSample &sample,
                   const TemporalGeneratedPaletteState &prepared, float hue) {
  const auto &state = prepared.current;
  return Color4(
      PaletteRgbHueLut::sample(prepared.lab, palette_value(sample, state), hue),
      palette_opacity(sample, state));
}

/** @brief Finishes a cached RGB phase through its flash-resident entry point.
 * @param sample Field coordinates and coverage.
 * @param prepared Cached interval base and palette controls.
 * @param hue Accumulated hue in turns.
 * @return Final straight-alpha color. */
HS_HOT_FLASH_MEMBER inline Color4
finish_cached_phase(const FieldSample &sample,
                    const TemporalGeneratedPaletteState &prepared, float hue) {
  return cached_phase_color(sample, prepared, hue);
}

/** @brief Finishes a cached RGB palette through one prepared phase.
 * @tparam Phase Prepared modulation phase type.
 * @tparam SinePhase Whether directional-sine phase sampling is supported.
 * @param sample Shared field coordinates.
 * @param prepared Prepared palette endpoints.
 * @param modulation Ready column phase.
 * @param output Receives the final straight-alpha color.
 * @return Whether this phase and RGB cache admit direct finishing. */
template <typename Phase, bool SinePhase = true>
HS_HOT_FLASH_MEMBER inline bool
shade_cached_phase(const FieldSample &sample,
                   const TemporalGeneratedPaletteState &prepared,
                   const Phase &modulation, Color4 &output)
  requires std::is_base_of_v<PhaseModulationState, Phase>
{
  if (cached_color_path(prepared, modulation, prepared.fraction,
                        modulation.ready) != CachedColorPath::RGB_PHASE)
    return false;
  PhaseScratch scratch;
  const float hue =
      modulation_active(modulation)
          ? accumulate_hue(0.0f,
                           Color::sample_phase_scalar<Phase, SinePhase>(
                               sample, modulation, scratch),
                           modulation.mechanic)
          : 0.0f;
  output = finish_cached_phase(sample, prepared, hue);
  return true;
}

/** @brief Accumulates one deferred positive lightness lift.
 * @param lift Prior lift.
 * @param scalar Source value.
 * @param mechanic Lift controls.
 * @return Combined unit-range lift. */
HS_HOT_INLINE inline float
accumulate_lift(float lift, float scalar,
                const ColorModulation::Params &mechanic) {
  const float lift_value = mechanic.amount * hs::clamp(scalar, 0.0f, 1.0f);
  return lift + (1.0f - lift) * lift_value;
}

/** @brief Updates deferred perceptual state from one scalar.
 * @param sample Carrier to update.
 * @param scalar Source value.
 * @param mechanic Ordered perceptual operation. */
HS_HOT_INLINE inline void
apply_perceptual_scalar(ColorSample &sample, float scalar,
                        const ColorModulation::Params &mechanic) {
  HS_AUDIT_CHECK(mechanic.kind != ColorModulation::Mechanic::OPACITY);
  if (mechanic.kind == ColorModulation::Mechanic::HUE_ROTATION)
    sample.pending_hue = accumulate_hue(sample.pending_hue, scalar, mechanic);
  else
    sample.pending_lift =
        accumulate_lift(sample.pending_lift, scalar, mechanic);
}

/** @brief Evaluates a bound phase source using its concrete prepared type.
 * @tparam Phase Type of the bound phase.
 * @tparam SinePhase Whether directional sine is supported.
 * @param field Spatial sample.
 * @param state Bound source phase of type Phase.
 * @param scratch Shared spatial lookups.
 * @return Signed scalar before color modulation. */
template <typename Phase, bool SinePhase>
HS_HOT_FLASH_MEMBER inline float
sample_bound_phase(const FieldSample &field, const PhaseModulationState &state,
                   PhaseScratch &scratch) {
  if constexpr (std::same_as<Phase, CubePhaseModulationState>) {
    if (state.source != Pattern::SourceKind::NOISE ||
        !static_cast<const CubePhaseModulationState &>(state).cube_current)
      return sample_bound_phase<PhaseModulationState, SinePhase>(field, state,
                                                                 scratch);
  }
  if constexpr (std::same_as<Phase, ChamberPhaseModulationState>)
    if (state.source != Pattern::SourceKind::NOISE)
      return sample_bound_phase<PhaseModulationState, SinePhase>(field, state,
                                                                 scratch);
  return sample_phase_scalar<Phase, SinePhase>(
      field, static_cast<const Phase &>(state), scratch);
}

/** @brief Applies a prepared column phase to a color carrier.
 * @tparam Phase Prepared modulation phase type.
 * @tparam SinePhase Whether directional-sine phase sampling is supported.
 * @param sample Carrier to update.
 * @param prepared Layer state.
 * @param scratch Cached analytic samples for this phase. */
template <typename Phase, bool SinePhase>
HS_HOT_FLASH_MEMBER inline void apply_phase_modulation(ColorSample &sample,
                                                       const Phase &prepared,
                                                       PhaseScratch &scratch)
  requires std::is_base_of_v<PhaseModulationState, Phase>
{
  if (!modulation_active(prepared))
    return;
  apply_perceptual_scalar(
      sample,
      sample_bound_phase<Phase, SinePhase>(sample.field, prepared, scratch),
      prepared.mechanic);
}

/** @brief Applies one modulation mechanic to a color carrier.
 * @param sample Carrier to update.
 * @param scalar Source value.
 * @param mechanic Ordered mechanic. */
HS_HOT_INLINE inline void
apply_mechanic_scalar(ColorSample &sample, float scalar,
                      const ColorModulation::Params &mechanic) {
  if (mechanic.kind == ColorModulation::Mechanic::HUE_ROTATION ||
      mechanic.kind == ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT)
    apply_perceptual_scalar(sample, scalar, mechanic);
  else {
    HS_CHECK(mechanic.kind == ColorModulation::Mechanic::OPACITY,
             "unknown modulation mechanic");
    sample.color = ColorModulation::Opacity::apply(sample.color, scalar,
                                                   {mechanic.amount});
  }
}
/** @brief Applies one ordered modulation layer.
 * @param sample Color and source coordinates, updated in place.
 * @param state Layer controls and endpoint sources.
 * @param previous Whether to sample the retained endpoint.
 * @tparam MixPath Whether the layer can blend path and spatial signals.
 */
template <bool MixPath = true>
HS_HOT_ITCM_MEMBER inline void
apply_modulation_in_place(ColorSample &sample, const ModulationState &state,
                          bool previous = false) {
  if (!modulation_active(state))
    return;
  float scalar;
  if (state.source == Pattern::SourceKind::DISPLACEMENT)
    scalar = sample.field.path_length;
  else {
    const auto *source = previous ? state.previous : state.current;
    HS_AUDIT_CHECK(source);
    scalar = Modulation::sample(
        *source, Pattern::SphereDirection::admit(sample.field.sphere));
  }
  apply_mechanic_scalar(sample,
                        state.scalar<MixPath>(scalar, sample.field.path_length),
                        state.mechanic);
}

/** @brief Applies one layer to retained and current color endpoints.
 * @tparam MixPath Whether the layer can blend path and spatial signals.
 * @param previous Retained color carrier.
 * @param current Current color carrier.
 * @param state Layer controls and immutable sources. */
template <bool MixPath = true>
HS_HOT_ITCM_MEMBER inline void
apply_modulation_pair_in_place(ColorSample &previous, ColorSample &current,
                               const ModulationState &state) {
  if (!modulation_active(state))
    return;
  if (state.source == Pattern::SourceKind::DISPLACEMENT ||
      previous.field.sphere != current.field.sphere) {
    apply_modulation_in_place<MixPath>(previous, state, true);
    apply_modulation_in_place<MixPath>(current, state, false);
    return;
  }
  HS_AUDIT_CHECK(state.previous && state.current);
  const auto direction = Pattern::SphereDirection::admit(current.field.sphere);
  const auto scalars =
      Modulation::sample_pair(*state.previous, *state.current, direction);
  apply_mechanic_scalar(
      previous,
      state.scalar<MixPath>(scalars.first, previous.field.path_length),
      state.mechanic);
  apply_mechanic_scalar(
      current, state.scalar<MixPath>(scalars.second, current.field.path_length),
      state.mechanic);
}

/** @brief Applies one ordered modulation layer to a copy.
 * @tparam MixPath Whether the layer can blend path and spatial signals.
 * @param input Color and source coordinates.
 * @param state Layer controls and endpoint sources.
 * @param previous Whether to sample the retained endpoint.
 * @return Modulated color carrier. */
template <bool MixPath = true>
inline ColorSample apply_modulation(const ColorSample &input,
                                    const ModulationState &state,
                                    bool previous = false) {
  ColorSample output = input;
  apply_modulation_in_place<MixPath>(output, state, previous);
  return output;
}

} // namespace Pullback::Color
