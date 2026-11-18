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
        count(sizeof...(Layers)), cached_base(cached) {
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
    const auto controls = palette_controls(field, palette.current);
    const auto previous =
        apply_palette_endpoint(field, palette.previous, controls.value,
                               controls.opacity, controls.gain);
    auto current =
        apply_palette_endpoint(field, palette.current, controls.value,
                               controls.opacity, controls.gain);
    blend_color_samples_in_place(previous, current, fraction);
    return shade_layers(current);
  }
};

} // namespace Pullback::Color
