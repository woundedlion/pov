/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file phase_stack.h
 * @brief Reuse of equal modulation sources at a shared column phase. */

#include "core/render/pullback/modulation_stack.h"
#include "core/render/pullback/palette_cache.h"
#include <tuple>

namespace Pullback::Color {

/** @brief One ordered layer's spatial sampling policy and column state.
 * @tparam Phase Prepared phase type.
 * @tparam Sine Whether this layer samples directional sine phases. */
template <typename Phase, bool Sine> struct PhaseLayer {
  using Prepared = Phase;                  ///< Typed column state.
  static constexpr bool SINE_PHASE = Sine; ///< Directional sine support.
  const Phase &state;                      ///< Prepared source and mechanic.
};

/** @brief Column-bound palette and ordered phase modulation sources. */
class PhaseModulationStack {
  struct Operation {
    const PhaseModulationState *state;
    float (*sample)(const FieldSample &, const PhaseModulationState &,
                    PhaseScratch &);
    uint8_t source;
  };
  using Apply = void (*)(ColorSample &, const PhaseModulationStack &);
  std::array<Operation, MAX_ORDERED_MODULATION_LAYERS> operations;
  const TemporalGeneratedPaletteState &palette;
  float fraction;
  Apply apply;
  uint8_t count;
  bool cached_base;
  uint16_t color_weight;
  float common_alpha = -1.0f;

  HS_HOT_FLASH_MEMBER ColorSample
  canonical_base(const FieldSample &field) const {
    const auto controls = palette_controls(field, palette.current);
    if (fraction <= 0.0f || fraction >= 1.0f)
      return apply_palette_endpoint_impl(
          field, fraction <= 0.0f ? palette.previous : palette.current,
          controls.value, controls.opacity, controls.gain);

    const auto &a = palette.previous;
    const auto &b = palette.current;
    const auto sample = [&](const BakedPalette &endpoint) {
      return common_alpha >= 0.0f
                 ? Color4(endpoint.get_color(controls.value), common_alpha)
                 : endpoint.get(controls.value);
    };
    const Color4 current = sample(*b.palette);
    const bool same_palette = a.palette == b.palette;
    const Color4 previous = same_palette ? current : sample(*a.palette);
    const bool same_rgb = previous.color == current.color;
    const auto *a_hue = a.hue_rotation.active ? a.hue_rotation.data : nullptr;
    const auto *b_hue = b.hue_rotation.active ? b.hue_rotation.data : nullptr;
    const bool retain_hue =
        same_rgb && a_hue == b_hue &&
        a.hue_rotation.gamut_mask == b.hue_rotation.gamut_mask;

    const Pixel rgb = same_rgb
                          ? current.color
                          : previous.color.lerp16(current.color, color_weight);
    const float alpha = current.alpha * controls.opacity;
    ColorSample result{
        Color4(rgb, same_palette || common_alpha >= 0.0f
                        ? alpha
                        : hs::lerp(previous.alpha * controls.opacity, alpha,
                                   fraction)),
        field,
        retain_hue ? b_hue : nullptr,
        controls.value,
        0.0f,
        controls.gain};
    result.hue_rotation_gamut =
        retain_hue ? b.hue_rotation.gamut_mask : nullptr;
    return result;
  }

  template <size_t I, typename Tuple, size_t... Earlier>
  void bind_source(const Tuple &layers, std::index_sequence<Earlier...>) {
    const auto &layer = std::get<I>(layers);
    using Layer = std::tuple_element_t<I, Tuple>;
    auto &op = operations[I];
    op = {&layer.state,
          sample_bound_phase<typename Layer::Prepared, Layer::SINE_PHASE>, I};
    if constexpr (I > 0) {
      if (!modulation_active(layer.state))
        return;
      const auto check = [&]<size_t J>() {
        const auto &prior = std::get<J>(layers);
        if (op.source == I && modulation_active(prior.state) &&
            same_modulation_phase<Layer::SINE_PHASE,
                                  std::tuple_element_t<J, Tuple>::SINE_PHASE>(
                layer.state, prior.state)) {
          op.source = operations[J].source;
          apply = &apply_repeated;
        }
      };
      (check.template operator()<Earlier>(), ...);
    }
  }

  template <typename... Layers>
  static void apply_distinct(ColorSample &sample,
                             const PhaseModulationStack &stack) {
    PhaseScratch scratch;
    [&]<size_t... I>(std::index_sequence<I...>) {
      (apply_phase_modulation<typename Layers::Prepared, Layers::SINE_PHASE>(
           sample,
           static_cast<const typename Layers::Prepared &>(
               *stack.operations[I].state),
           scratch),
       ...);
    }(std::index_sequence_for<Layers...>{});
  }

  HS_HOT_FLASH_MEMBER static void
  apply_repeated(ColorSample &sample, const PhaseModulationStack &stack) {
    PhaseScratch scratch;
    std::array<float, MAX_ORDERED_MODULATION_LAYERS> values;
    for (uint8_t i = 0; i < stack.count; ++i) {
      const auto &op = stack.operations[i];
      if (modulation_active(*op.state)) {
        if (op.source == i)
          values[i] = op.sample(sample.field, *op.state, scratch);
        apply_perceptual_scalar(sample, values[op.source], op.state->mechanic);
      }
    }
  }

  static bool active_byte_cube(const CubePhaseModulationState &layer) {
    return modulation_active(layer) &&
           layer.source == Pattern::SourceKind::NOISE && layer.cube_previous &&
           layer.cube_current && layer.cube_previous->data &&
           layer.cube_current->data && !layer.cube_previous->linear_data &&
           !layer.cube_current->linear_data;
  }

  HS_HOT_ITCM_MEMBER static void
  apply_byte_cubes(ColorSample &sample, const PhaseModulationStack &stack) {
    const auto tap =
        Pattern::prepare_sphere_lut_sample_unchecked(sample.field.sphere);
    std::array<float, MAX_ORDERED_MODULATION_LAYERS> values;
    for (uint8_t i = 0; i < stack.count; ++i) {
      const auto &op = stack.operations[i];
      if (op.source == i) {
        const auto &cube =
            static_cast<const CubePhaseModulationState &>(*op.state);
        const auto *a = cube.cube_previous->data;
        const auto *b = cube.cube_current->data;
        const float previous = Pattern::bounded_scalar(tap.sample_data(a));
        values[i] = tap.same_byte_corners(a, b)
                        ? previous
                        : hs::lerp(previous,
                                   Pattern::bounded_scalar(tap.sample_data(b)),
                                   cube.cube_fraction);
      }
      apply_perceptual_scalar(sample, values[op.source], op.state->mechanic);
    }
  }

  HS_HOT_INLINE Color4 shade_layers(ColorSample &sample) const {
    apply(sample, *this);
    return finish_color_sample(sample);
  }

public:
  /** @brief Binds phase sources and palette cache admission in authored order.
   * @tparam Layers Typed phase layers.
   * @param base Palette endpoints.
   * @param weight Current column weight.
   * @param cached Whether the interval base cache admits this column.
   * @param inputs Prepared phase layer views. */
  template <typename... Layers>
  PhaseModulationStack(const TemporalGeneratedPaletteState &base, float weight,
                       bool cached, Layers... inputs)
      : palette(base), fraction(weight), apply(&apply_distinct<Layers...>),
        count(sizeof...(Layers)), cached_base(cached),
        color_weight(frac_to_q16(weight)) {
    if (!cached && weight > 0.0f && weight < 1.0f &&
        base.previous.palette != base.current.palette) {
      const auto a = base.previous.palette->constant_alpha();
      const auto b = base.current.palette->constant_alpha();
      if (a && b == a)
        common_alpha = *a;
    }
    static_assert(sizeof...(Layers) <= MAX_ORDERED_MODULATION_LAYERS);
    const auto layers = std::tuple{inputs...};
    [&]<size_t... I>(std::index_sequence<I...>) {
      (bind_source<I>(layers, std::make_index_sequence<I>{}), ...);
    }(std::index_sequence_for<Layers...>{});
    if constexpr ((std::same_as<typename Layers::Prepared,
                                CubePhaseModulationState> &&
                   ...))
      if ((active_byte_cube(inputs.state) && ...))
        apply = &apply_byte_cubes;
  }

  /** @brief Shades a field using its shared palette base and bound phase sources.
   * @param field Pixel coordinates and coverage.
   * @return Final straight-alpha color. */
  HS_HOT_ITCM_MEMBER Color4 shade(const FieldSample &field) const {
    if (cached_base) {
      ColorSample sample;
      GeneratedPaletteCore<true, true, false>::apply_shared_base_ready(
          field, palette, sample);
      return shade_layers(sample);
    }
    auto sample = canonical_base(field);
    return shade_layers(sample);
  }
};

} // namespace Pullback::Color
