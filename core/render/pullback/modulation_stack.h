/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/pullback/color_state.h"
#include <array>
#include <cstring>
#include <span>

/** @file modulation_stack.h
 * @brief Shared sampling of ordered endpoint color modulations. */

namespace Pullback::Color {

/// Maximum number of ordered layers in the bounded column kernel.
inline constexpr size_t MAX_ORDERED_MODULATION_LAYERS = 3;

namespace StackDetail {
/** @brief Identifies the effective interpolation data of two borrowed cubes. */
inline bool same_cube(const Pattern::SphereLutView &a,
                      const Pattern::SphereLutView &b) {
  return a.linear_data == b.linear_data && (a.linear_data || a.data == b.data);
}

/** @brief Obtains the additional state of an exact phase representation. */
template <typename T, typename Phase>
constexpr const T *phase_as(const Phase &phase) {
  if constexpr (std::is_same_v<T, Phase>)
    return &phase;
  else
    return nullptr;
}
} // namespace StackDetail

/** @brief Checks whether two prepared column phases produce identical scalars.
 * @tparam LeftSine Whether the first sampler supports directional sine.
 * @tparam RightSine Whether the second sampler supports directional sine.
 * @param a First prepared phase.
 * @param b Second prepared phase.
 * @return Whether one sample can serve both mechanics for this column. */
template <bool LeftSine, bool RightSine, typename Left, typename Right>
HS_SIZE_FLASH_MEMBER bool same_modulation_phase(const Left &a, const Right &b) {
  if (!a.ready || !b.ready || a.source != b.source)
    return false;
  switch (a.source) {
  case Pattern::SourceKind::DISPLACEMENT:
    return true;
  case Pattern::SourceKind::DIRECTIONAL_SINE:
    return LeftSine && RightSine &&
           std::memcmp(&a.grid.a, &b.grid.a, sizeof(a.grid.a)) == 0;
  case Pattern::SourceKind::SPHERE_SMOOTH_GRID:
    return std::memcmp(&a.grid, &b.grid, sizeof(a.grid)) == 0;
  case Pattern::SourceKind::NOISE:
    break;
  default:
    return false;
  }
  const auto *cube_a = StackDetail::phase_as<CubePhaseModulationState>(a);
  const auto *cube_b = StackDetail::phase_as<CubePhaseModulationState>(b);
  if ((cube_a && cube_a->cube_current) || (cube_b && cube_b->cube_current))
    return cube_a && cube_b && cube_a->cube_previous && cube_b->cube_previous &&
           cube_a->cube_current && cube_b->cube_current &&
           cube_a->cube_fraction == cube_b->cube_fraction &&
           StackDetail::same_cube(*cube_a->cube_previous,
                                  *cube_b->cube_previous) &&
           StackDetail::same_cube(*cube_a->cube_current, *cube_b->cube_current);
  const auto *patch_a = StackDetail::phase_as<ChamberPhaseModulationState>(a);
  const auto *patch_b = StackDetail::phase_as<ChamberPhaseModulationState>(b);
  const bool cached_a =
      patch_a && patch_a->noise_previous && patch_a->noise_current;
  const bool cached_b =
      patch_b && patch_b->noise_previous && patch_b->noise_current;
  if (cached_a != cached_b ||
      (cached_a &&
       (patch_a->noise_previous != patch_b->noise_previous ||
        patch_a->noise_current != patch_b->noise_current ||
        patch_a->noise_fraction != patch_b->noise_fraction ||
        std::memcmp(&patch_a->noise_to_chamber, &patch_b->noise_to_chamber,
                    sizeof(patch_a->noise_to_chamber)) != 0)))
    return false;
  const bool lazy_a = patch_a && patch_a->lazy_noise_offset;
  const bool lazy_b = patch_b && patch_b->lazy_noise_offset;
  if (lazy_a != lazy_b ||
      (lazy_a ? patch_a->noise_phase != patch_b->noise_phase
              : std::memcmp(&a.noise_offset, &b.noise_offset,
                            sizeof(a.noise_offset)) != 0))
    return false;
  return a.noise && b.noise &&
         (a.noise == b.noise ||
          (a.noise->scale == b.noise->scale &&
           a.noise->basis == b.noise->basis &&
           a.noise->variant == b.noise->variant &&
           a.noise->noise.same_noise_single_3d(b.noise->noise)));
}

/** @brief Prepared endpoint sources and authored-order mechanic instructions. */
class OrderedModulationKernel {
  static constexpr uint8_t SOURCE_CAPACITY = 2 * MAX_ORDERED_MODULATION_LAYERS;
  static constexpr uint8_t PATH = SOURCE_CAPACITY;
  enum HueLiftMaskSlot : uint8_t {
    HUE,
    LIFT,
    PREVIOUS_LIFT,
    MASK,
    PREVIOUS_MASK
  };
  enum class SourceKind : uint8_t { DIRECT, BYTE_CUBE, FLOAT_CUBE };
  struct Source {
    const Modulation::Prepared *value;
    SourceKind kind;
    uint8_t peer;
  };
  struct Operation {
    ColorModulation::Params mechanic;
    uint8_t current;
    uint8_t previous;
  };

  const GeneratedPaletteState &previous_palette;
  float fraction;
  const GeneratedPaletteState &base;
  bool same_palette;
  bool blend;
  bool retained;
  bool identity_controls;
  bool same_hue_table;
  uint16_t rgb_weight;
  bool valid = false;
  bool uses_cube = false;
  bool hue_lift_mask = false;
  uint8_t source_count = 0;
  uint8_t op_count = 0;
  std::array<Source, SOURCE_CAPACITY> sources;
  std::array<Operation, MAX_ORDERED_MODULATION_LAYERS> ops{};

  static bool identity(const GeneratedPaletteState &state) {
    return state.mapping.exact ==
               static_cast<uint8_t>(PaletteMapping::LINEAR) &&
           state.mapping_frequency == 1.0f && state.mapping_offset == 0.0f &&
           state.brightness_envelope == BrightnessEnvelope::NONE &&
           state.opacity_low == 1.0f && state.opacity_high == 1.0f;
  }

  static SourceKind kind(const Modulation::Prepared &endpoint) {
    if (endpoint.evaluate != Modulation::sample_cube || !endpoint.cube.active)
      return SourceKind::DIRECT;
    return endpoint.cube.linear_data ? SourceKind::FLOAT_CUBE
                                     : SourceKind::BYTE_CUBE;
  }
  static bool can_sample(const Modulation::Prepared *source) {
    return source && source->evaluate &&
           (source->evaluate != Modulation::sample_cube ||
            (source->cube.active &&
             (source->cube.data || source->cube.linear_data)));
  }

  HS_SIZE_FLASH_MEMBER uint8_t bind_source(const Modulation::Prepared &value,
                                           uint8_t peer = PATH) {
    const auto type = kind(value);
    for (uint8_t i = 0; i < source_count; ++i) {
      const auto &other = *sources[i].value;
      if (&value == &other)
        return i;
      if (sources[i].kind != type)
        continue;
      if (type == SourceKind::BYTE_CUBE && value.cube.data == other.cube.data)
        return i;
      if (type == SourceKind::FLOAT_CUBE &&
          value.cube.linear_data == other.cube.linear_data)
        return i;
      if (type == SourceKind::DIRECT &&
          Modulation::same_direct_sample(value, other))
        return i;
    }
    const uint8_t i = source_count++;
    HS_AUDIT_CHECK(i < SOURCE_CAPACITY);
    sources[i] = {&value, type, peer};
    uses_cube |= type != SourceKind::DIRECT;
    return i;
  }

  bool byte_pair(const Operation &op) const {
    return op.current != PATH && op.previous != PATH &&
           sources[op.current].kind == SourceKind::BYTE_CUBE &&
           sources[op.previous].kind == SourceKind::BYTE_CUBE;
  }

  HS_SIZE_FLASH_MEMBER bool bind_hue_lift_mask() {
    if (op_count != 3)
      return false;
    const auto &[hue, lift, mask] = ops;
    using Mechanic = ColorModulation::Mechanic;
    const bool matches =
        hue.mechanic.kind == Mechanic::HUE_ROTATION &&
        lift.mechanic.kind == Mechanic::POSITIVE_LIGHTNESS_LIFT &&
        mask.mechanic.kind == Mechanic::OPACITY && hue.current != PATH &&
        hue.current == hue.previous &&
        sources[hue.current].kind == SourceKind::DIRECT && byte_pair(lift) &&
        byte_pair(mask) && lift.current != mask.current &&
        lift.current != mask.previous && lift.previous != mask.current &&
        lift.previous != mask.previous;
    if (!matches)
      return false;
    const std::array endpoints = {
        sources[hue.current].value, sources[lift.current].value,
        sources[lift.previous].value, sources[mask.current].value,
        sources[mask.previous].value};
    for (size_t i = 0; i < endpoints.size(); ++i)
      sources[i].value = endpoints[i];
    return true;
  }

  static HS_HOT_INLINE float sample_bytes(const Pattern::SphereLutSample &tap,
                                          const Pattern::SphereLutView &cube) {
    return Pattern::bounded_scalar(tap.sample_data(cube.data));
  }

  static HS_HOT_INLINE float paired_bytes(const Pattern::SphereLutSample &tap,
                                          const Pattern::SphereLutView &cube,
                                          const Pattern::SphereLutView &peer,
                                          float peer_value) {
    const auto *a = cube.data + tap.offset;
    const auto *b = peer.data + tap.offset;
    constexpr int STRIDE = Pattern::SphereLutView::FACE_STEPS;
    if (std::memcmp(a, b, 2) == 0 &&
        std::memcmp(a + STRIDE, b + STRIDE, 2) == 0)
      return peer_value;
    return sample_bytes(tap, cube);
  }

  static HS_HOT_INLINE float opacity(float alpha, float scalar, float amount) {
    return ColorModulation::Opacity::apply(Color4(Pixel(0, 0, 0), alpha),
                                           scalar, {amount})
        .alpha;
  }

  /** @brief Blends palette RGB first and retains separate endpoint alpha. */
  HS_HOT_FLASH_MEMBER ColorSample palette_base(const FieldSample &field,
                                               float &previous_alpha) const {
    auto result = identity_controls
                      ? apply_palette_endpoint(field, base, field.value,
                                               field.coverage, 1.0f)
                      : apply_generated_palette(field, base);
    previous_alpha = result.color.alpha;
    if (blend) {
      bool keep_hue_table = same_hue_table;
      if (!same_palette) {
        const auto old =
            previous_palette.palette->get(result.palette_coordinate);
        previous_alpha =
            old.alpha *
            (identity_controls ? field.coverage : palette_opacity(field, base));
        keep_hue_table &= old.color == result.color.color;
        result.color.color = old.color.lerp16(result.color.color, rgb_weight);
      }
      if (!keep_hue_table) {
        result.hue_rotation = nullptr;
        result.hue_rotation_gamut = nullptr;
      }
    }
    return result;
  }

  HS_HOT_FLASH_MEMBER Color4
  shade_hue_lift_mask(const FieldSample &field) const {
    const auto &[hue, lift, mask] = ops;
    float old_alpha;
    auto result = palette_base(field, old_alpha);
    const auto direction = Pattern::SphereDirection::admit(field.sphere);
    result.pending_hue = accumulate_hue(
        0.0f, Modulation::sample(*sources[HUE].value, direction), hue.mechanic);
    const auto tap = Pattern::prepare_sphere_lut_sample(direction);
    const auto &lift_cube = sources[LIFT].value->cube;
    const auto &mask_cube = sources[MASK].value->cube;
    const float l = sample_bytes(tap, lift_cube);
    const float m = sample_bytes(tap, mask_cube);
    result.pending_lift = accumulate_lift(0.0f, l, lift.mechanic);
    result.color.alpha = opacity(result.color.alpha, m, mask.mechanic.amount);
    if (blend) {
      const float old_l =
          paired_bytes(tap, sources[PREVIOUS_LIFT].value->cube, lift_cube, l);
      const float old_m =
          paired_bytes(tap, sources[PREVIOUS_MASK].value->cube, mask_cube, m);
      result.pending_lift =
          hs::lerp(accumulate_lift(0.0f, old_l, lift.mechanic),
                   result.pending_lift, fraction);
      result.color.alpha =
          hs::lerp(opacity(old_alpha, old_m, mask.mechanic.amount),
                   result.color.alpha, fraction);
    }
    return finish_color_sample(result);
  }

  HS_HOT_FLASH_MEMBER Color4 shade_ordered(const FieldSample &field) const {
    float previous_alpha;
    auto result = palette_base(field, previous_alpha);
    std::array<float, SOURCE_CAPACITY + 1> values;
    values[PATH] = field.path_length;
    if (source_count) {
      const auto direction = Pattern::SphereDirection::admit(field.sphere);
      Pattern::SphereLutSample tap{};
      if (uses_cube)
        tap = Pattern::prepare_sphere_lut_sample(direction);
      for (uint8_t i = 0; i < source_count; ++i) {
        const auto &source = sources[i];
        const auto &value = *source.value;
        switch (source.kind) {
        case SourceKind::DIRECT:
          values[i] = Modulation::sample(value, direction);
          break;
        case SourceKind::BYTE_CUBE:
          if (source.peer != PATH &&
              sources[source.peer].kind == SourceKind::BYTE_CUBE)
            values[i] =
                paired_bytes(tap, value.cube, sources[source.peer].value->cube,
                             values[source.peer]);
          else
            values[i] = sample_bytes(tap, value.cube);
          break;
        case SourceKind::FLOAT_CUBE:
          values[i] =
              Pattern::bounded_scalar(tap.sample_data(value.cube.linear_data));
          break;
        }
      }
    }
    float previous_hue = 0.0f;
    float previous_lift = 0.0f;
    for (uint8_t i = 0; i < op_count; ++i) {
      const auto &op = ops[i];
      const float scalar = values[op.current];
      switch (op.mechanic.kind) {
      case ColorModulation::Mechanic::HUE_ROTATION:
        result.pending_hue =
            accumulate_hue(result.pending_hue, scalar, op.mechanic);
        if (blend)
          previous_hue =
              accumulate_hue(previous_hue, values[op.previous], op.mechanic);
        break;
      case ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT:
        result.pending_lift =
            accumulate_lift(result.pending_lift, scalar, op.mechanic);
        if (blend)
          previous_lift =
              accumulate_lift(previous_lift, values[op.previous], op.mechanic);
        break;
      case ColorModulation::Mechanic::OPACITY:
        result.color.alpha =
            opacity(result.color.alpha, scalar, op.mechanic.amount);
        if (blend)
          previous_alpha =
              opacity(previous_alpha, values[op.previous], op.mechanic.amount);
        break;
      }
    }
    if (blend) {
      result.pending_hue =
          blend_unit_hue(previous_hue, result.pending_hue, fraction);
      result.pending_lift =
          hs::lerp(previous_lift, result.pending_lift, fraction);
      result.color.alpha =
          hs::lerp(previous_alpha, result.color.alpha, fraction);
    }
    return finish_color_sample(result);
  }

protected:
  /** @brief Binds active layers in authored order and interns distinct sources.
   * @param current Current palette.
   * @param previous Retained palette.
   * @param t Current endpoint weight.
   * @param same Whether both palette endpoints have equal entries.
   * @param layers Consecutive canonical modulation stages. */
  HS_SIZE_FLASH_MEMBER
  OrderedModulationKernel(const GeneratedPaletteState &current,
                          const GeneratedPaletteState &previous, float t,
                          bool same,
                          std::span<const ModulationState *const> layers)
      : previous_palette(previous), fraction(t),
        base(t <= 0.0f ? previous : current), same_palette(same),
        blend(t > 0.0f && t < 1.0f), retained(t <= 0.0f),
        identity_controls(identity(base)),
        same_hue_table((previous.hue_rotation.active
                            ? previous.hue_rotation.data
                            : nullptr) == (current.hue_rotation.active
                                               ? current.hue_rotation.data
                                               : nullptr) &&
                       previous.hue_rotation.gamut_mask ==
                           current.hue_rotation.gamut_mask),
        rgb_weight(blend && !same ? frac_to_q16(t) : 0) {
    if (layers.size() > MAX_ORDERED_MODULATION_LAYERS || !base.palette ||
        (blend && !previous.palette))
      return;
    for (const auto *layer : layers) {
      if (!layer)
        return;
      if (!modulation_active(*layer))
        continue;
      if (layer->mechanic.kind > ColorModulation::Mechanic::OPACITY ||
          layer->source > Pattern::SourceKind::DISPLACEMENT)
        return;
      uint8_t current = PATH, previous = PATH;
      if (layer->source != Pattern::SourceKind::DISPLACEMENT) {
        const auto *selected = retained ? layer->previous : layer->current;
        if (!can_sample(selected) || (blend && !can_sample(layer->previous)))
          return;
        current = bind_source(*selected);
        previous = blend && layer->previous != selected
                       ? bind_source(*layer->previous, current)
                       : current;
      }
      ops[op_count++] = {layer->mechanic, current, previous};
    }
    valid = true;
  }

  template <size_t N> void bind_specialization() {
    if constexpr (N == 3)
      hue_lift_mask = valid && bind_hue_lift_mask();
  }

  template <size_t N>
  HS_HOT_INLINE Color4 shade_bound(const FieldSample &field) const {
    HS_CHECK(ready(), "invalid ordered modulation column");
    if constexpr (N == 3)
      if (hue_lift_mask)
        return shade_hue_lift_mask(field);
    return shade_ordered(field);
  }

public:
  /** @brief Whether all active layers were bound.
   * @return False when the canonical pipeline must handle the column. */
  bool ready() const { return valid; }
};

/** @brief Column binding for a typed sequence of modulation stages.
 * @tparam N Number of consecutive stages. */
template <size_t N> struct ModulationStack : OrderedModulationKernel {
  /** @brief Binds the ordered modulation stages to a generated palette.
   * @param palette Prepared temporal palette.
   * @param same_palette Whether the palette entries match.
   * @param layers Prepared layer endpoints. */
  ModulationStack(const TemporalGeneratedPaletteState &palette,
                  bool same_palette,
                  std::array<const ModulationState *, N> layers)
      : OrderedModulationKernel(palette.current, palette.previous,
                                palette.fraction, same_palette, layers) {
    bind_specialization<N>();
  }

  /** @brief Binds consecutive layers on one non-temporal palette.
   * @param palette Current palette.
   * @param layers Prepared layer endpoints. */
  ModulationStack(const GeneratedPaletteState &palette,
                  std::array<const ModulationState *, N> layers)
      : OrderedModulationKernel(palette, palette, 1.0f, true, layers) {
    bind_specialization<N>();
  }

  /** @brief Evaluates the bound column at one field sample.
   * @param field Shared source direction, coverage and displacement.
   * @return Finished straight-alpha color. */
  HS_HOT_INLINE Color4 shade(const FieldSample &field) const {
    return shade_bound<N>(field);
  }
};

} // namespace Pullback::Color
