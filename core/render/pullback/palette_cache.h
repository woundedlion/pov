/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Generated-palette cache selection and reconstruction. */

#include "core/render/pullback/color_state.h"
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

/// Linear-RGB distance from a gamut gate requiring the canonical palette base.
inline constexpr float PALETTE_LAB_GAMUT_GUARD_WIDTH =
    PaletteLabLutBase::GAMUT_GUARD_WIDTH;

/** @brief Resolves lightness and hue into rotated perceptual coordinates.
 * @param lab Base perceptual color.
 * @param hue Hue rotation in turns.
 * @param lift Lightness lift toward white.
 * @return Perceptual color before gamut mapping. */
HS_HOT_INLINE inline OKLab rotated_palette_lab(OKLab lab, float hue,
                                               float lift) {
  lab.L += (1.0f - lab.L) * lift;
  float cosine;
  float sine;
  hue_cos_sin(hue, cosine, sine);
  return rotate_chroma(lab, cosine, sine);
}

/** @brief Tests raw RGB proximity to a discontinuous gamut gate.
 * @param raw Hue-rotated linear RGB before gamut mapping.
 * @return Whether the canonical palette base is required. */
HS_HOT_INLINE inline bool near_palette_gamut_gate(const LinRGB &raw) {
  return PaletteLabLutBase::near_gamut_gate(raw);
}

/** @brief Restores a cached base whose raw RGB can cross the gamut gate.
 * @tparam Restore Canonical base callback.
 * @param sample Deferred color carrier, updated on fallback.
 * @param restore Replaces the carrier's RGB base and hue-table identity.
 * @return Whether the canonical callback ran. */
template <typename Restore>
HS_HOT_INLINE inline bool
restore_cached_palette_base_near_gamut_gate(ColorSample &sample,
                                            Restore &&restore) {
  if (!sample.has_cached_lab || sample.cached_rgb ||
      sample.brightness_gain == 0.0f)
    return false;
  if (!near_palette_gamut_gate(oklab_to_linear_rgb(rotated_palette_lab(
          sample.cached_lab, sample.pending_hue, sample.pending_lift))))
    return false;
  restore();
  sample.has_cached_lab = false;
  return true;
}

/** @brief Finishes perceptual color outside the canonical gamut gate.
 * @param lab Cached perceptual base.
 * @param hue Deferred hue rotation in turns.
 * @param lift Deferred lightness lift.
 * @param gain Final RGB brightness gain.
 * @param alpha Straight alpha.
 * @param result Receives the finished color when the gamut gate is clear.
 * @return Whether the cached base supplied the finished color. */
HS_HOT_ITCM_MEMBER inline bool try_finish_palette_lab(OKLab lab, float hue,
                                                      float lift, float gain,
                                                      float alpha,
                                                      Color4 &result) {
  if (gain == 0.0f) {
    result = Color4(Pixel(0, 0, 0), alpha);
    return true;
  }
  const auto lab_value = rotated_palette_lab(lab, hue, lift);
  const auto raw = oklab_to_linear_rgb(lab_value);
  if (near_palette_gamut_gate(raw))
    return false;
  auto output = gamut_map_lut(raw, lab_value);
  if (gain != 1.0f) {
    apply_brightness_gain(output, gain);
  }
  result = Color4(linrgb_to_pixel(output), alpha);
  return true;
}

/** @brief Finishes perceptual color with canonical gamut-gate fallback.
 * @tparam Canonical Canonical finished-color callback.
 * @param lab Cached perceptual base.
 * @param hue Deferred hue rotation in turns.
 * @param lift Deferred lightness lift.
 * @param gain Final RGB brightness gain.
 * @param alpha Straight alpha.
 * @param canonical Produces canonical RGB-first finished color on fallback.
 * @return Final straight-alpha color. */
template <typename Canonical>
HS_HOT_INLINE inline Color4 finish_palette_lab(OKLab lab, float hue, float lift,
                                               float gain, float alpha,
                                               Canonical &&canonical) {
  Color4 result;
  if (try_finish_palette_lab(lab, hue, lift, gain, alpha, result))
    return result;
  return canonical();
}

/** @brief Finishes a guarded cached carrier using one raw RGB evaluation.
 * @tparam Restore Canonical base callback.
 * @param sample Deferred color carrier, updated on fallback.
 * @param restore Replaces the carrier's RGB base and hue-table identity.
 * @return Final straight-alpha color. */
template <typename Restore>
HS_HOT_INLINE inline Color4 finish_fused_color(ColorSample &sample,
                                               Restore &&restore) {
  if (!sample.has_cached_lab || sample.cached_rgb)
    return finish_color_sample(sample);
  return finish_palette_lab(sample.cached_lab, sample.pending_hue,
                            sample.pending_lift, sample.brightness_gain,
                            sample.color.alpha, [&] {
                              restore();
                              sample.has_cached_lab = false;
                              return finish_color_sample(sample);
                            });
}

/// Qualified hue-rotation and noise-LUT modulation bounds.
inline constexpr std::array<ApproximationMetric, 3> MODULATION_LUT_METRICS{{
    {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MAXIMUM,
     7000.0f, "channel code"},
    {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MEAN, 256.0f,
     "channel code"},
    {ApproximationDomain::FRAMEBUFFER, ApproximationAggregation::MAXIMUM,
     3066.0f, "channel code"},
}};
static_assert(HueRotationGamutMask::MAX_CHANNEL_SPAN ==
              2 * MODULATION_LUT_METRICS[0].limit);

/** @brief Samples canonical modulation endpoints at shared coordinates.
 * @param sample Shared field coordinates.
 * @param modulation Prepared modulation endpoints.
 * @return Previous and current endpoint scalars. */
template <bool MixPath = true>
HS_HOT_ITCM_MEMBER inline std::pair<float, float>
sample_fused_scalars(const FieldSample &sample,
                     const ModulationState &modulation) {
  if (modulation.source == Pattern::SourceKind::DISPLACEMENT)
    return {sample.path_length, sample.path_length};
  HS_AUDIT_CHECK(modulation.previous && modulation.current);
  const auto direction = Pattern::SphereDirection::admit(sample.sphere);
  if (modulation.paired_cube_ready) {
    const auto tap = Pattern::prepare_sphere_lut_sample(direction);
    return {Pattern::bounded_scalar(tap.sample(modulation.previous->cube)),
            Pattern::bounded_scalar(tap.sample(modulation.current->cube))};
  }
  const auto pair = Modulation::sample_pair(*modulation.previous,
                                            *modulation.current, direction);
  if constexpr (MixPath)
    return {modulation.scalar<MixPath>(pair.first, sample.path_length),
            modulation.scalar<MixPath>(pair.second, sample.path_length)};
  else
    return pair;
}

/** @brief Restores the canonical RGB-first palette base.
 * @param result Deferred carrier to restore.
 * @param prepared Prepared palette endpoints.
 * @param fraction Current endpoint weight. */
HS_HOT_ITCM_MEMBER inline void
restore_palette_base(ColorSample &result,
                     const TemporalGeneratedPaletteState &prepared,
                     float fraction) {
  const Pixel a =
      prepared.previous.palette->get_color(result.palette_coordinate);
  const Pixel b =
      prepared.current.palette->get_color(result.palette_coordinate);
  const Pixel *const a_table = prepared.previous.hue_rotation.active
                                   ? prepared.previous.hue_rotation.data
                                   : nullptr;
  const Pixel *const b_table = prepared.current.hue_rotation.active
                                   ? prepared.current.hue_rotation.data
                                   : nullptr;
  result.color.color = a.lerp16(b, frac_to_q16(fraction));
  result.hue_rotation = a == b && a_table == b_table ? b_table : nullptr;
  result.hue_rotation_gamut =
      result.hue_rotation ? prepared.current.hue_rotation.gamut_mask : nullptr;
}

/** @brief Palette sampling shared by endpoint and cache configuration.
 * @tparam TEMPORAL Whether previous palettes are retained.
 * @tparam LAB_CACHE Whether an RGB-first base cache is available.
 * @tparam FUSED_TEMPORAL_COLOR Whether endpoint modulation may be fused. */
template <bool TEMPORAL, bool LAB_CACHE, bool FUSED_TEMPORAL_COLOR>
struct GeneratedPaletteCore {
  using Prepared =
      std::conditional_t<TEMPORAL, TemporalGeneratedPaletteState,
                         GeneratedPaletteState>; ///< Endpoint state.
  /** @brief Produces one color endpoint.
   * @param sample Field coordinates.
   * @param prepared Prepared endpoints.
   * @param previous Selects the retained endpoint.
   * @return Color carrier. */
  static ColorSample apply_endpoint(const FieldSample &sample,
                                    const Prepared &prepared, bool previous) {
    if constexpr (TEMPORAL)
      return apply_generated_palette(sample, previous ? prepared.previous
                                                      : prepared.current);
    else {
      (void)previous;
      return apply_generated_palette(sample, prepared);
    }
  }
  /** @brief Samples both palettes with their shared mapping and gain controls.
   * @param sample Shared field coordinates.
   * @param prepared Prepared endpoints.
   * @return Retained endpoint followed by current endpoint. */
  static std::pair<ColorSample, ColorSample>
  apply_pair(const FieldSample &sample, const Prepared &prepared)
    requires TEMPORAL
  {
    const auto &state = prepared.current;
    const auto controls = palette_controls(sample, state);
    if constexpr (LAB_CACHE) {
      const auto endpoint = [&](const GeneratedPaletteState &endpoint_state)
                                HS_HOT_INLINE {
                                  return apply_cache_fallback_endpoint(
                                      sample, endpoint_state, controls.value,
                                      controls.opacity, controls.gain);
                                };
      return {endpoint(prepared.previous), endpoint(state)};
    } else {
      return {apply_palette_endpoint(sample, prepared.previous, controls.value,
                                     controls.opacity, controls.gain),
              apply_palette_endpoint(sample, state, controls.value,
                                     controls.opacity, controls.gain)};
    }
  }
  /** @brief Samples a canonical endpoint for a cache-enabled fallback.
   * @param sample Shared field coordinates.
   * @param state Endpoint resources.
   * @param value Prepared palette coordinate.
   * @param opacity Shared coverage and opacity.
   * @param gain Shared final brightness.
   * @return Canonical deferred endpoint. */
  HS_HOT_FLASH_MEMBER static ColorSample
  apply_cache_fallback_endpoint(const FieldSample &sample,
                                const GeneratedPaletteState &state, float value,
                                float opacity, float gain)
    requires LAB_CACHE
  {
    return apply_palette_endpoint_impl(sample, state, value, opacity, gain);
  }
  /** @brief Checks the cache view and temporal fraction for a shared base.
   * @param prepared Prepared endpoints and cache view.
   * @return Whether the cache supplies this column's palette base. */
  static bool shared_base_ready(const Prepared &prepared) {
    if constexpr (TEMPORAL && LAB_CACHE)
      return prepared.lab.low &&
             !(prepared.fraction <= 0.0f || prepared.fraction > 1.0f);
    else
      return false;
  }
  /** @brief Tries one cached base without fetching either palette.
   * @param sample Shared field coordinates.
   * @param prepared Prepared endpoint controls and cache view.
   * @param output Receives the cached color carrier.
   * @return Whether the cache supplied the base. */
  HS_HOT_INLINE static bool apply_shared_base(const FieldSample &sample,
                                              const Prepared &prepared,
                                              ColorSample &output)
    requires(TEMPORAL && LAB_CACHE)
  {
    if (!shared_base_ready(prepared))
      return false;
    apply_shared_base_ready(sample, prepared, output);
    return true;
  }

  /** @brief Produces a cached base after column eligibility was checked.
   * @param sample Shared field coordinates.
   * @param prepared Prepared eligible cache view.
   * @param output Receives the cached color carrier.
   * @pre The base view is active and its fraction is inside (0, 1]. */
  HS_HOT_INLINE static void apply_shared_base_ready(const FieldSample &sample,
                                                    const Prepared &prepared,
                                                    ColorSample &output)
    requires(TEMPORAL && LAB_CACHE)
  {
    const auto &state = prepared.current;
    const float value = palette_value(sample, state);
    output.color = Color4(Pixel(0, 0, 0), palette_opacity(sample, state));
    output.field = sample;
    output.hue_rotation = nullptr;
    output.palette_coordinate = value;
    output.pending_hue = 0.0f;
    output.brightness_gain = palette_gain(sample, state);
    output.pending_lift = 0.0f;
    output.cached_rgb =
        prepared.lab.rgb_low && prepared.lab.rgb_high ? &prepared.lab : nullptr;
    output.cached_lab = output.cached_rgb
                            ? OKLab{}
                            : PaletteLabLut::sample(prepared.lab, value);
    output.has_cached_lab = true;
  }

private:
  /** @brief Resolves one canonical deferred hue endpoint.
   * @param scalar Source scalar.
   * @param amount Hue turns per source unit.
   * @return Deferred hue rotation in turns. */
  HS_HOT_INLINE static float fused_endpoint_hue(float scalar, float amount) {
    return accumulate_hue(0.0f, scalar, {.amount = amount});
  }

  /** @brief Blends wrapped endpoint hues along their shortest path.
   * @param previous Retained endpoint hue.
   * @param current Current endpoint hue.
   * @param fraction Current endpoint weight.
   * @return Blended hue in turns. */
  HS_HOT_INLINE static float blend_fused_hue(float previous, float current,
                                             float fraction) {
    return blend_unit_hue(previous, current, fraction);
  }

  /** @brief Resolves the canonical shortest-path endpoint hue blend.
   * @param scalars Previous and current source scalars.
   * @param amount Hue turns per source unit.
   * @param fraction Current endpoint weight.
   * @return Deferred hue rotation in turns. */
  HS_HOT_INLINE static float
  fused_pending_hue(const std::pair<float, float> &scalars, float amount,
                    float fraction) {
    const float previous = fused_endpoint_hue(scalars.first, amount);
    const float current = fused_endpoint_hue(scalars.second, amount);
    return blend_fused_hue(previous, current, fraction);
  }

  /** @brief Blends palette and modulation endpoints into one color carrier.
   * @param result Fresh color carrier to populate.
   * @param sample Shared field coordinates.
   * @param prepared Prepared palette endpoints.
   * @param modulation Prepared modulation endpoints.
   * @param fraction Current endpoint weight. */
  template <bool MixPath>
  HS_HOT_INLINE static void
  prepare_fused_modulation(ColorSample &result, const FieldSample &sample,
                           const Prepared &prepared,
                           const ModulationState &modulation, float fraction)
    requires(TEMPORAL && FUSED_TEMPORAL_COLOR)
  {
    if (fraction <= 0.0f || fraction >= 1.0f) {
      result = apply_endpoint(sample, prepared, fraction <= 0.0f);
      apply_modulation_in_place<MixPath>(result, modulation, fraction <= 0.0f);
      return;
    }
    const bool cached = [&] {
      if constexpr (LAB_CACHE)
        return apply_shared_base(sample, prepared, result);
      else
        return false;
    }();
    Color4 previous;
    Color4 current;
    const Pixel *previous_table;
    const Pixel *current_table;
    if (cached) {
      previous = current = result.color;
      previous_table = current_table = nullptr;
    } else {
      const auto &state = prepared.current;
      const float value = palette_value(sample, state);
      const float opacity = palette_opacity(sample, state);
      previous = prepared.previous.palette->get(value);
      current = state.palette->get(value);
      previous.alpha *= opacity;
      current.alpha *= opacity;
      previous_table = prepared.previous.hue_rotation.active
                           ? prepared.previous.hue_rotation.data
                           : nullptr;
      current_table =
          state.hue_rotation.active ? state.hue_rotation.data : nullptr;
      result.field = sample;
      result.hue_rotation = current_table;
      result.hue_rotation_gamut = prepared.current.hue_rotation.gamut_mask;
      result.palette_coordinate = value;
      result.brightness_gain = palette_gain(sample, state);
    }
    float previous_hue = 0.0f;
    float current_hue = 0.0f;
    float previous_lift = 0.0f;
    float current_lift = 0.0f;
    if (modulation_active(modulation)) {
      const auto scalars = sample_fused_scalars<MixPath>(sample, modulation);
      const auto apply = [&](Color4 &color, float &hue, float &lift,
                             float scalar) {
        if (modulation.mechanic.kind == ColorModulation::Mechanic::HUE_ROTATION)
          hue = accumulate_hue(hue, scalar, modulation.mechanic);
        else if (modulation.mechanic.kind ==
                 ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT) {
          lift = accumulate_lift(lift, scalar, modulation.mechanic);
        } else
          color = ColorModulation::apply(color, scalar, modulation.mechanic);
      };
      apply(previous, previous_hue, previous_lift, scalars.first);
      apply(current, current_hue, current_lift, scalars.second);
    }
    const bool same_base =
        previous.color == current.color && previous_table == current_table;
    result.color =
        Color4(previous.color.lerp16(current.color, frac_to_q16(fraction)),
               hs::lerp(previous.alpha, current.alpha, fraction));
    result.pending_hue = blend_fused_hue(previous_hue, current_hue, fraction);
    result.pending_lift = hs::lerp(previous_lift, current_lift, fraction);
    if (!same_base)
      result.hue_rotation = nullptr;
  }

  /** @brief Restores the canonical RGB-first palette base.
   * @param result Deferred carrier to restore.
   * @param prepared Prepared palette endpoints.
   * @param fraction Current endpoint weight. */
  HS_HOT_INLINE static void restore_fused_palette_base(ColorSample &result,
                                                       const Prepared &prepared,
                                                       float fraction) {
    Color::restore_palette_base(result, prepared, fraction);
  }

public:
  /** @brief Blends palette and modulation endpoints into one guarded carrier.
   * @param sample Shared field coordinates.
   * @param prepared Prepared palette endpoints.
   * @param modulation Prepared modulation endpoints.
   * @param fraction Current endpoint weight.
   * @return Deferred color after canonical endpoint blending. */
  template <bool MixPath = true>
  HS_HOT_ITCM_MEMBER static ColorSample
  apply_fused_modulation(const FieldSample &sample, const Prepared &prepared,
                         const ModulationState &modulation, float fraction)
    requires(TEMPORAL && FUSED_TEMPORAL_COLOR)
  {
    ColorSample result;
    prepare_fused_modulation<MixPath>(result, sample, prepared, modulation,
                                      fraction);
    if constexpr (LAB_CACHE)
      restore_cached_palette_base_near_gamut_gate(result, [&] {
        restore_fused_palette_base(result, prepared, fraction);
      });
    return result;
  }

  /** @brief Shades fused endpoints without returning a deferred carrier.
   * @param sample Shared field coordinates.
   * @param prepared Prepared palette endpoints.
   * @param modulation Prepared modulation endpoints.
   * @param fraction Current endpoint weight.
   * @return Final straight-alpha color. */
  template <bool MixPath = true>
  HS_HOT_FLASH_MEMBER static Color4
  shade_fused_modulation(const FieldSample &sample, const Prepared &prepared,
                         const ModulationState &modulation, float fraction)
    requires(TEMPORAL && FUSED_TEMPORAL_COLOR)
  {
    if constexpr (LAB_CACHE)
      switch (cached_color_path(prepared, modulation, fraction, false)) {
      case CachedColorPath::LAB_ENDPOINT: {
        const auto &state = prepared.current;
        const float value = palette_value(sample, state);
        const float opacity = palette_opacity(sample, state);
        const float hue = [&] {
          if (!modulation_active(modulation))
            return 0.0f;
          const float scalar = [&] {
            if (modulation.source == Pattern::SourceKind::DISPLACEMENT)
              return sample.path_length;
            HS_AUDIT_CHECK(modulation.current);
            return Modulation::sample(
                *modulation.current,
                Pattern::SphereDirection::admit(sample.sphere));
          }();
          return fused_endpoint_hue(
              modulation.scalar<MixPath>(scalar, sample.path_length),
              modulation.mechanic.amount);
        }();
        return finish_palette_lab(
            PaletteLabLut::sample(prepared.lab, value), hue, 0.0f, 1.0f,
            opacity, [&] {
              ColorSample canonical{
                  Color4(state.palette->get_color(value), opacity), {}};
              canonical.palette_coordinate = value;
              canonical.pending_hue = hue;
              canonical.hue_rotation =
                  state.hue_rotation.active ? state.hue_rotation.data : nullptr;
              canonical.hue_rotation_gamut = state.hue_rotation.gamut_mask;
              return finish_color_sample(canonical);
            });
      }
      case CachedColorPath::LAB_HUE_INTERVAL: {
        const auto scalars = sample_fused_scalars<MixPath>(sample, modulation);
        const auto &state = prepared.current;
        const float value = palette_value(sample, state);
        const float opacity = palette_opacity(sample, state);
        const float hue =
            fused_pending_hue(scalars, modulation.mechanic.amount, fraction);
        return finish_palette_lab(
            PaletteLabLut::sample(prepared.lab, value), hue, 0.0f, 1.0f,
            opacity, [&] {
              ColorSample canonical{Color4(Pixel(0, 0, 0), opacity), {}};
              canonical.palette_coordinate = value;
              canonical.pending_hue = hue;
              restore_fused_palette_base(canonical, prepared, fraction);
              return finish_color_sample(canonical);
            });
      }
      case CachedColorPath::CANONICAL:
      case CachedColorPath::RGB_PHASE:
        break;
      }
    ColorSample result;
    prepare_fused_modulation<MixPath>(result, sample, prepared, modulation,
                                      fraction);
    if constexpr (LAB_CACHE)
      return finish_fused_color(result, [&] {
        restore_fused_palette_base(result, prepared, fraction);
      });
    else
      return finish_color_sample(result);
  }

  /** @brief Produces the current color endpoint.
   * @param sample Field coordinates.
   * @param prepared Prepared endpoints.
   * @return Color carrier. */
  static ColorSample apply(const FieldSample &sample,
                           const Prepared &prepared) {
    return apply_endpoint(sample, prepared, false);
  }
};

} // namespace Pullback::Color
