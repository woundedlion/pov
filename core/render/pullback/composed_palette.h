/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>
#include <variant>
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
#if HS_ENABLE_PARAM_GUI_BRIDGE
  static constexpr bool PALETTE_SOURCE_FADE = AUTHORED_PALETTE;

  struct SourceFade {
    AuthoredPaletteResources target;
    BakedPaletteStorage mixed;
    PaletteCycler target_cycle;
    uint32_t target_hue = 0;
    float target_chroma = -1.0f;
    const ColorFamily *from = nullptr;
    const ColorFamily *to = nullptr;
    float progress = 0.0f;
    uint32_t generation = 0;
    bool initialized = false;
  };
#endif
#if HS_ENABLE_TEST_HOOKS
  friend struct hs_test::ComposedFrameWhiteBox;
#endif

public:
  /** @brief Persistent byte budget for all lazy palette fade resources.
   * @return Upper bound including the alternate cycle and both baked tables. */
  static constexpr size_t source_fade_arena_bytes() {
#if HS_ENABLE_PARAM_GUI_BRIDGE
    if constexpr (PALETTE_SOURCE_FADE)
      return sizeof(SourceFade) + alignof(SourceFade) +
             PaletteCycler::generated_arena_bytes() +
             2 * BakedPalette::required_arena_bytes();
#endif
    return 0;
  }

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
#if HS_ENABLE_PARAM_GUI_BRIDGE
    clear_source_fade();
#endif
  }

  /** @brief Captures the mapping endpoints of an arming crossfade.
   * @param target Crossfade's destination color controls. */
  HS_O3_FN HS_HOT_INLINE void arm(const ColorFamily &target) {
    mapping_from = mapping;
    mapping_to = Color::PaletteMappingWeights::single(target.palette_mapping);
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /** @brief Drops any armed palette-source endpoints. */
  HS_O3_FN HS_HOT_INLINE void clear_source_fade() {
    if constexpr (PALETTE_SOURCE_FADE)
      if (source_fade)
        source_fade->from = nullptr;
  }

  /** @brief Retains the controls of the active preset edge.
   * @param from Outgoing palette controls.
   * @param target Incoming palette controls. */
  HS_O3_FN HS_HOT_INLINE void arm_source(const ColorFamily &from,
                                         const ColorFamily &target) {
    if constexpr (PALETTE_SOURCE_FADE) {
      clear_source_fade();
      if (different_selections(from, target)) {
        if (!source_fade)
          source_fade = persistent_arena.make<SourceFade>();
        source_fade->from = &from;
        source_fade->to = &target;
        source_fade->progress = 0.0f;
        if (cycling(from) && cycling(target))
          arm_target_cycle();
      }
    }
  }

#endif

  /** @brief Interpolates the mapping weights of an in-flight crossfade.
   * @param progress Eased transition progress, 0 to 1. */
  HS_O3_FN HS_HOT_INLINE void blend(float progress) {
    mapping =
        Color::PaletteMappingWeights::lerp(mapping_from, mapping_to, progress);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    if constexpr (PALETTE_SOURCE_FADE)
      if (source_fade && source_fade->from)
        source_fade->progress = progress;
#endif
  }

  /** @brief Steps the palette mapping wobble phase. */
  HS_O3_FN HS_HOT_INLINE void advance() {
    oscillation_phase =
        math::wrap_t(oscillation_phase + color->phase_oscillation_speed);
  }

  /** @brief Follows the chroma and authored cycling harmony controls. */
  HS_FLASH_MEMBER void update_chroma() {
#if HS_ENABLE_PARAM_GUI_BRIDGE
    const float chroma = cycling_chroma();
#else
    const float chroma = color->palette_chroma;
#endif
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
#if HS_ENABLE_PARAM_GUI_BRIDGE
      if (!cycling(cycling_controls())) {
#else
      if (!cycling(*color)) {
#endif
        cycler.advance_without_display();
        return;
      }
    }
    cycler.step();
#if HS_ENABLE_PARAM_GUI_BRIDGE
    if constexpr (PALETTE_SOURCE_FADE)
      if (cycling_harmony_fade()) {
        auto &fade = *source_fade;
        if (fade.target_chroma != chroma) {
          fade.target_chroma = chroma;
          fade.target_cycle.defer_generated_chroma(chroma);
        }
        fade.target_cycle.step();
      }
#endif
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
  HS_SIZE_FLASH_MEMBER void prepare(FrameState &frame, Storage &storage) {
#if HS_ENABLE_PARAM_GUI_BRIDGE
    const auto [palette_ptr, palette_generation] = prepared_palette(storage);
    const BakedPalette &palette = *palette_ptr;
#else
    const BakedPalette &palette = selected_palette(storage);
    uint32_t palette_generation = cycler.bake_generation() << 1;
    if constexpr (AUTHORED_PALETTE)
      if (!cycling(*color))
        palette_generation =
            (storage.authored_palette.bake_generation() << 1) | 1U;
#endif
    if constexpr (HUE_ROTATION_LUT)
      storage.refresh_hue_rotation_lut(palette, palette_generation);
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
      if (!cycling(*color))
        return storage.authored_palette.prepare(persistent_arena, *color,
                                                color->palette_chroma);
    }
    return cycler.palette();
  }

private:
  /** @brief Identifies the cycling palette source. */
  static bool cycling(const ColorFamily &controls) {
    return controls.palette_type ==
           static_cast<uint8_t>(AuthoredPaletteType::CYCLING_GENERATED);
  }

  /** @brief Resolves a palette mode to its generated harmony. */
  HS_O3_FN HS_HOT_INLINE static PaletteHarmony harmony_for_mode(uint8_t mode) {
    constexpr auto &HARMONIES = Color::PALETTE_MODE_HARMONIES;
    HS_CHECK(mode < std::size(HARMONIES), "invalid cycling palette mode");
    return HARMONIES[mode];
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /** @brief Checks whether a discrete authored control changes the source. */
  static bool different_selections(const ColorFamily &from,
                                   const ColorFamily &to) {
    if (from.palette_type != to.palette_type)
      return true;
    using Type = AuthoredPaletteType;
    switch (static_cast<Type>(from.palette_type)) {
    case Type::CYCLING_GENERATED:
      return from.palette_mode != to.palette_mode;
    case Type::GENERATED:
      return from.palette_harmony != to.palette_harmony ||
             from.palette_lightness_curve != to.palette_lightness_curve;
    case Type::PROCEDURAL:
      return from.palette_name != to.palette_name;
    case Type::MUTATING:
      return from.palette_name != to.palette_name ||
             from.palette_name_b != to.palette_name_b;
    default:
      return false;
    }
  }

  /** @brief Whether both retained endpoints use different cycling harmonies. */
  bool cycling_harmony_fade() const {
    if constexpr (PALETTE_SOURCE_FADE)
      return source_fade && source_fade->from && cycling(*source_fade->from) &&
             cycling(*source_fade->to);
    return false;
  }

  /** @brief Reads the controls of the cycling endpoint, if it contributes. */
  const ColorFamily &cycling_controls() const {
    if constexpr (PALETTE_SOURCE_FADE)
      if (source_fade && source_fade->from) {
        if (cycling(*source_fade->from))
          return *source_fade->from;
        if (cycling(*source_fade->to))
          return *source_fade->to;
      }
    return *color;
  }

  /** @brief Chroma shared by the generated cycles. */
  float cycling_chroma() const {
    return cycling_harmony_fade() ? color->palette_chroma
                                  : cycling_controls().palette_chroma;
  }

  /** @brief Starts or resynchronizes the target harmony at the current clock. */
  HS_FLASH_MEMBER void arm_target_cycle() {
    auto &fade = *source_fade;
    const bool initialized = fade.target_chroma >= 0.0f;
    fade.target_chroma = cycling_chroma();
    if (!initialized)
      fade.target_cycle.init_generated(
          persistent_arena, {.next = next_target_palette,
                             .context = this,
                             .dwell_frames = GeneratedPaletteBank::DWELL_FRAMES,
                             .fade_frames = GeneratedPaletteBank::FADE_FRAMES,
                             .easing = math::ease_in_out_sin,
                             .display_interval_frames = 2});
    fade.target_hue = hue;
    GeneratedPaletteBank::restore_cycler(
        fade.target_cycle, cycler.generated_clock(), fade.target_hue,
        harmony_for_mode(fade.to->palette_mode), fade.target_chroma);
  }

  /** @brief Generates a palette for the incoming harmony. */
  static void next_target_palette(void *context, uint32_t sequence,
                                  GenerativePalette &out) {
    auto &fade = *static_cast<ComposedPalette *>(context)->source_fade;
    GeneratedPaletteBank::next_palette(fade.target_hue, sequence,
                                       harmony_for_mode(fade.to->palette_mode),
                                       fade.target_chroma, out);
  }

  /** @brief Resolves one authored endpoint against the shared cycle clock. */
  const BakedPalette &resolve(AuthoredPaletteResources &storage,
                              const ColorFamily &controls) {
    if (cycling(controls))
      return cycler.palette();
    return storage.prepare(persistent_arena, controls, controls.palette_chroma);
  }

  /** @brief Returns the frame's palette and hue-table cache key. */
  template <typename Storage>
  auto prepared_palette(Storage &storage)
      -> std::pair<const BakedPalette *, uint32_t> {
    constexpr uint32_t TAG_BITS = PALETTE_SOURCE_FADE ? 2 : 1;
    if constexpr (PALETTE_SOURCE_FADE) {
      auto *fade = source_fade;
      if (fade && fade->from && fade->progress > 0.0f) {
        const auto &from = resolve(storage.authored_palette, *fade->from);
        const auto &to = cycling_harmony_fade()
                             ? fade->target_cycle.palette()
                             : resolve(fade->target, *fade->to);
        if (fade->initialized)
          fade->mixed.rebake_crossfade(from, to, fade->progress);
        else {
          fade->mixed.bake_blend(persistent_arena, from, to, fade->progress);
          fade->initialized = true;
        }
        return {&fade->mixed.view(), (++fade->generation << TAG_BITS) | 2U};
      }
    }
    if constexpr (AUTHORED_PALETTE)
      if (!cycling(*color)) {
        const auto &palette = storage.authored_palette.prepare(
            persistent_arena, *color, color->palette_chroma);
        return {&palette,
                (storage.authored_palette.bake_generation() << TAG_BITS) | 1U};
      }
    return {&cycler.palette(), cycler.bake_generation() << TAG_BITS};
  }

#endif

  /** @brief Current generated-cycle harmony. @return Selected hue relationship. */
  HS_FLASH_MEMBER PaletteHarmony cycling_harmony() const {
    if constexpr (AUTHORED_PALETTE) {
#if HS_ENABLE_PARAM_GUI_BRIDGE
      return harmony_for_mode(cycling_controls().palette_mode);
#else
      return harmony_for_mode(color->palette_mode);
#endif
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
#if HS_ENABLE_PARAM_GUI_BRIDGE
    GeneratedPaletteBank::next_palette(palette.hue, sequence,
                                       palette.cycling_harmony(),
                                       palette.cycling_chroma(), out);
#else
    GeneratedPaletteBank::next_palette(palette.hue, sequence,
                                       palette.cycling_harmony(),
                                       palette.color->palette_chroma, out);
#endif
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
#if HS_ENABLE_PARAM_GUI_BRIDGE
#if defined(_MSC_VER)
  [[msvc::no_unique_address]]
#else
  [[no_unique_address]]
#endif
  std::conditional_t<PALETTE_SOURCE_FADE, SourceFade *, std::monostate>
      source_fade{};
#endif
};

} // namespace Pullback
