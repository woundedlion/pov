/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstdint>
#include <span>
#include "core/color/authored_palette.h"
#include "core/color/palette_cycler.h"
#include "core/memory.h"
#include "core/render/pullback/composed_descriptors.h"

#if HS_ENABLE_TEST_HOOKS
namespace hs_test {
struct ComposedFrameWhiteBox;
}
#endif

/** @file composed_palette.h
 * @brief `ComposedPalette`, the palette lifecycle of a composed effect.
 */

namespace Pullback {

/**
 * @brief Generated palette cycle, palette mapping weights, mapping wobble
 *        and the per-frame palette and hue-table resolution of one composed
 *        effect.
 * @details Reads the effect's live color controls. Arena-owned palette
 * storage (`authored_palette` and the `OptionalHueRotationLut` members) stays
 * in the runtime State and is passed to the calls that read it.
 * @tparam Traits The effect's `ComposedTraits`.
 */
template <typename Traits> class ComposedPalette {
  using ColorFamily = typename Traits::ColorFamily;
  using FrameState = typename Traits::FrameState;
  using SpecT = typename Traits::Spec;
  static constexpr bool AUTHORED_PALETTE = Traits::AUTHORED_PALETTE;
  static constexpr bool HUE_ROTATION_LUT = Traits::HUE_ROTATION_LUT;
#if HS_ENABLE_TEST_HOOKS
  friend struct hs_test::ComposedFrameWhiteBox;
#endif

public:
  /** @brief Binds the effect's live color controls.
   * @param color Color controls the palette reads for its lifetime. */
  HS_HOT_INLINE explicit ComposedPalette(const ColorFamily &color)
      : color(&color),
        mapping(Color::PaletteMappingWeights::single(color.palette_mapping)) {}

  /** @brief Starts the generated palette cycle in `persistent_arena`. */
  HS_O3_FN HS_HOT_INLINE void init() {
    cycler.init_generated(persistent_arena,
                          {.next = next_palette,
                           .context = this,
                           .dwell_frames = GeneratedPaletteBank::DWELL_FRAMES,
                           .fade_frames = GeneratedPaletteBank::FADE_FRAMES,
                           .easing = math::ease_in_out_sin,
                           .display_interval_frames = 2});
  }

  /** @brief Re-derives the mapping weights from adopted controls.
   * @param target Adopted color controls. */
  HS_O3_FN HS_HOT_INLINE void adopt(const ColorFamily &target) {
    mapping = Color::PaletteMappingWeights::single(target.palette_mapping);
  }

  /** @brief Captures the mapping endpoints of an arming crossfade.
   * @param target Crossfade's destination color controls. */
  HS_O3_FN HS_HOT_INLINE void arm(const ColorFamily &target) {
    mapping_from = mapping;
    mapping_to = Color::PaletteMappingWeights::single(target.palette_mapping);
  }

  /** @brief Interpolates the mapping weights of an in-flight crossfade.
   * @param progress Eased transition progress, 0 to 1. */
  HS_O3_FN HS_HOT_INLINE void blend(float progress) {
    mapping =
        Color::PaletteMappingWeights::lerp(mapping_from, mapping_to, progress);
  }

  /** @brief Steps the palette mapping wobble phase. */
  HS_O3_FN HS_HOT_INLINE void advance() {
    oscillation_phase =
        math::wrap_t(oscillation_phase + color->phase_oscillation_speed);
  }

  /** @brief Follows the chroma and authored cycling harmony controls. */
  HS_FLASH_MEMBER void update_chroma() {
    const float chroma = color->palette_chroma;
    if constexpr (AUTHORED_PALETTE) {
      const auto harmony = cycling_harmony();
      if (this->harmony != harmony) {
        GeneratedPaletteBank::restore_cycler(cycler, cycler.generated_clock(),
                                             hue, harmony, chroma);
        this->harmony = harmony;
      }
    }
    if (this->chroma == chroma)
      return;
    this->chroma = chroma;
    cycler.defer_generated_chroma(this->chroma);
  }

  /** @brief Follows the color controls and steps the generated cycle. */
  HS_O3_FN HS_HOT_INLINE void step() {
    update_chroma();
    if constexpr (AUTHORED_PALETTE) {
      if (color->palette_type !=
          static_cast<uint8_t>(AuthoredPaletteType::CYCLING_GENERATED)) {
        cycler.advance_without_display();
        return;
      }
    }
    cycler.step();
  }

  /**
   * @brief Resolves this draw's palette, rebakes a stale hue table and
   *        writes the frame's palette fields.
   * @tparam Storage Runtime State holding `authored_palette` and the
   *         `OptionalHueRotationLut` members.
   * @param frame Receives the palette, hue table, mapping and wobble phase.
   * @param storage Arena-owned palette storage.
   */
  template <typename Storage>
  HS_O3_FN HS_HOT_INLINE void prepare(FrameState &frame, Storage &storage) {
    const BakedPalette &palette = selected_palette(storage);
    uint32_t palette_generation = cycler.bake_generation() << 1;
    if constexpr (AUTHORED_PALETTE)
      if (color->palette_type !=
          static_cast<uint8_t>(AuthoredPaletteType::CYCLING_GENERATED))
        palette_generation =
            (storage.authored_palette.bake_generation() << 1) | 1U;
    if constexpr (HUE_ROTATION_LUT)
      if (storage.hue_rotation_lut_bake != palette_generation) {
        Color::prepare_hue_rotation_lut(
            std::span<Pixel, Color::HueRotationLutView::SIZE>(
                storage.hue_rotation_lut),
            palette, &storage.hue_rotation_gamut);
        storage.hue_rotation_lut_bake = palette_generation;
      }
    frame.palette = &palette;
    if constexpr (HUE_ROTATION_LUT) {
      frame.hue_rotation_lut = storage.hue_rotation_lut.data();
      frame.hue_rotation_gamut = &storage.hue_rotation_gamut;
    } else {
      frame.hue_rotation_lut = nullptr;
      frame.hue_rotation_gamut = nullptr;
    }
    frame.palette_mapping = mapping;
    frame.palette_oscillation_phase = oscillation_phase;
  }

  /** @brief Resolves the current authored or cycling source.
   * @tparam Storage Runtime State holding `authored_palette`.
   * @param storage Arena-owned palette storage.
   * @return Palette table. */
  template <typename Storage>
  const BakedPalette &selected_palette(Storage &storage) {
    if constexpr (AUTHORED_PALETTE) {
      if (color->palette_type !=
          static_cast<uint8_t>(AuthoredPaletteType::CYCLING_GENERATED))
        return storage.authored_palette.prepare(persistent_arena, *color,
                                                color->palette_chroma);
    }
    return cycler.palette();
  }

private:
  /** @brief Current generated-cycle harmony. @return Selected hue relationship. */
  HS_FLASH_MEMBER PaletteHarmony cycling_harmony() const {
    if constexpr (AUTHORED_PALETTE) {
      constexpr auto &HARMONIES = Color::PALETTE_MODE_HARMONIES;
      const auto mode = color->palette_mode;
      HS_CHECK(mode < std::size(HARMONIES), "invalid cycling palette mode");
      return HARMONIES[mode];
    } else
      return SpecT::HARMONY;
  }

  /**
   * @brief Generator the palette cycler calls for each palette in the cycle.
   * @details Advances the hue on every palette after the first, so the cycle
   * opens on `hue`, which starts at 0 and moves only here.
   * @param context The ComposedPalette, as registered with the cycler.
   * @param sequence Zero-based index of the palette being generated.
   * @param out Receives the recipe to bake.
   */
  static void next_palette(void *context, uint32_t sequence,
                           GenerativePalette &out) {
    ComposedPalette &palette = *static_cast<ComposedPalette *>(context);
    GeneratedPaletteBank::next_palette(palette.hue, sequence,
                                       palette.cycling_harmony(),
                                       palette.color->palette_chroma, out);
  }

  const ColorFamily *color;                  ///< Live color controls.
  Color::PaletteMappingWeights mapping;      ///< Current mapping weights.
  Color::PaletteMappingWeights mapping_from; ///< Crossfade origin weights.
  Color::PaletteMappingWeights mapping_to;   ///< Crossfade target weights.
  float oscillation_phase = 0.0f;            ///< Mapping wobble phase.
  float chroma = -1.0f;                      ///< Chroma the cycle generates at.
  uint32_t hue = 0;                          ///< Generated-cycle hue clock.
  PaletteHarmony harmony = SpecT::HARMONY;   ///< Generated-cycle harmony.
  PaletteCycler cycler;                      ///< Generated palette cycle.
};

} // namespace Pullback
