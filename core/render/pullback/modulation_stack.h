/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/pullback/color_state.h"
#include <array>
#include <bit>
#include <cstring>
#include <limits>
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
    union {
      const Modulation::Prepared *value;
      const int8_t *bytes;
    };
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
  bool simple_palette;
  uint16_t rgb_weight;
  float constant_alpha = -1.0f;
  bool uses_cube = false;
  bool hue_lift_mask = false;
  Modulation::ImmutableTable::TapEvaluator hue_tap = nullptr;
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

  /** @brief Keeps unit directions unchanged when binary32 normalization divides by one. */
  static HS_HOT_INLINE Pattern::SphereDirection
  sample_direction(const math::Vector &v) {
    const float squared = v.x * v.x + v.y * v.y + v.z * v.z;
    // sqrt rounds to one at 1 and its immediate binary32 successor.
    if (squared >= 1.0f &&
        squared <= 1.0f + std::numeric_limits<float>::epsilon())
      return Pattern::SphereDirection::from_unit(v);
    if (squared < 1.0f &&
        squared >= 1.0f - std::numeric_limits<float>::epsilon()) {
      // Dividing a normal binary32 by the predecessor of one increments its
      // magnitude by exactly one representable step. Subnormals round unchanged.
      const auto divide = [](float component) HS_HOT_INLINE {
        const auto bits = std::bit_cast<uint32_t>(component);
        return std::bit_cast<float>(bits + ((bits & 0x7f800000u) != 0));
      };
      return Pattern::SphereDirection::from_unit(
          {divide(v.x), divide(v.y), divide(v.z)});
    }
    return Pattern::SphereDirection::admit(v);
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
    sources[i] = {{&value}, type, peer};
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
    const auto *direct_hue = sources[hue.current].value;
    if (direct_hue->immutable_table &&
        direct_hue->evaluate == direct_hue->immutable_table->evaluator)
      hue_tap = direct_hue->immutable_table->sample_tap;
    const std::array tables = {sources[lift.current].value->cube.data,
                               sources[lift.previous].value->cube.data,
                               sources[mask.current].value->cube.data,
                               sources[mask.previous].value->cube.data};
    sources[HUE].value = direct_hue;
    for (size_t i = 0; i < tables.size(); ++i)
      sources[LIFT + i].bytes = tables[i];
    if (simple_palette)
      constant_alpha = base.palette->constant_alpha().value_or(-1.0f);
    return true;
  }

  static HS_HOT_INLINE float sample_bytes(const Pattern::SphereLutSample &tap,
                                          const int8_t *cube) {
    return Pattern::bounded_scalar(tap.sample_data(cube));
  }

  static HS_HOT_INLINE float paired_bytes(const Pattern::SphereLutSample &tap,
                                          const int8_t *cube,
                                          const int8_t *peer,
                                          float peer_value) {
    if (tap.same_byte_corners(cube, peer))
      return peer_value;
    return sample_bytes(tap, cube);
  }

  /** @brief Samples the positive lobe of a signed-byte cube. */
  static HS_HOT_INLINE float positive_bytes(const Pattern::SphereLutSample &tap,
                                            const int8_t *cube) {
    const auto *row = cube + tap.offset;
    uint16_t low, high;
    std::memcpy(&low, row, sizeof(low));
    std::memcpy(&high, row + Pattern::SphereLutView::FACE_STEPS, sizeof(high));
    if ((low & high & 0x8080u) == 0x8080u)
      return 0.0f;
    return hs::clamp(tap.sample_data(cube), 0.0f, 1.0f);
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

  HS_HOT_INLINE float mask_alpha(float current, float previous, float scalar,
                                 float old_scalar, float amount) const {
    const float alpha = opacity(current, scalar, amount);
    return blend ? hs::lerp(opacity(previous, old_scalar, amount), alpha,
                            fraction)
                 : alpha;
  }

  HS_HOT_FLASH_MEMBER Color4 finish_blended_palette(const FieldSample &field,
                                                    float hue, float lift,
                                                    float mask,
                                                    float old_mask) const {
    float old_alpha;
    auto sample = palette_base(field, old_alpha);
    sample.pending_hue = hue;
    sample.pending_lift = lift;
    sample.color.alpha = mask_alpha(sample.color.alpha, old_alpha, mask,
                                    old_mask, ops[2].mechanic.amount);
    return finish_color_sample(sample);
  }

  HS_HOT_ITCM_MEMBER Color4
  shade_hue_lift_mask(const FieldSample &field) const {
    const auto &[hue, lift, mask] = ops;
    const auto direction = sample_direction(field.sphere);
    const auto tap = Pattern::prepare_sphere_lut_sample(direction);
    const auto &hue_source = *sources[HUE].value;
    const float value = hue_tap
                            ? Pattern::bounded_scalar(hue_tap(hue_source, tap))
                            : Modulation::sample(hue_source, direction);
    const float rotation = math::wrap_t(hue.mechanic.amount * value);
    const auto *lift_cube = sources[LIFT].bytes;
    const auto *mask_cube = sources[MASK].bytes;
    const float l = positive_bytes(tap, lift_cube);
    const float m = sample_bytes(tap, mask_cube);
    float lightness = lift.mechanic.amount * l;
    float old_mask = m;
    if (blend) {
      const auto *old_cube = sources[PREVIOUS_LIFT].bytes;
      const float old_l =
          old_cube == lift_cube ? l : positive_bytes(tap, old_cube);
      const auto *old_mask_cube = sources[PREVIOUS_MASK].bytes;
      old_mask =
          old_mask_cube == mask_cube ? m : sample_bytes(tap, old_mask_cube);
      lightness = hs::lerp(lift.mechanic.amount * old_l, lightness, fraction);
    }
    if (!simple_palette)
      return finish_blended_palette(field, rotation, lightness, m, old_mask);
    const float alpha =
        (constant_alpha >= 0.0f ? constant_alpha
                                : base.palette->get_alpha(field.value)) *
        field.coverage;
    const float result_alpha =
        mask_alpha(alpha, alpha, m, old_mask, mask.mechanic.amount);
    if (lightness == 0.0f && rotation != 0.0f && base.hue_rotation.active) {
      Pixel rotated;
      if (sample_guarded_hue_lut(base.hue_rotation.data,
                                 base.hue_rotation.gamut_mask, field.value,
                                 rotation, rotated))
        return Color4(rotated, result_alpha);
    }
    const Pixel color = base.palette->get_color(field.value);
    if (lightness == 0.0f && rotation == 0.0f)
      return Color4(color, result_alpha);
    return Color4(
        finish_hue_lift_pixel(color, nullptr, lightness, rotation, 1.0f),
        result_alpha);
  }

  HS_HOT_FLASH_MEMBER Color4 shade_ordered(const FieldSample &field) const {
    float previous_alpha;
    auto result = palette_base(field, previous_alpha);
    std::array<float, SOURCE_CAPACITY + 1> values;
    values[PATH] = field.path_length;
    if (source_count) {
      const auto direction = sample_direction(field.sphere);
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
            values[i] = paired_bytes(tap, value.cube.data,
                                     sources[source.peer].value->cube.data,
                                     values[source.peer]);
          else
            values[i] = sample_bytes(tap, value.cube.data);
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

  HS_HOT_INLINE bool
  bind_layers(std::span<const ModulationState *const> layers) {
    if (layers.size() > MAX_ORDERED_MODULATION_LAYERS || !base.palette ||
        (blend && !previous_palette.palette))
      return false;
    for (const auto *layer : layers) {
      if (!layer)
        return false;
      if (!modulation_active(*layer))
        continue;
      if (layer->mechanic.kind > ColorModulation::Mechanic::OPACITY ||
          layer->source > Pattern::SourceKind::DISPLACEMENT)
        return false;
      uint8_t current = PATH, previous = PATH;
      if (layer->source != Pattern::SourceKind::DISPLACEMENT) {
        const auto *selected = retained ? layer->previous : layer->current;
        if (!can_sample(selected) || (blend && !can_sample(layer->previous)))
          return false;
        current = bind_source(*selected);
        previous = blend && layer->previous != selected
                       ? bind_source(*layer->previous, current)
                       : current;
      }
      ops[op_count++] = {layer->mechanic, current, previous};
    }
    return true;
  }

protected:
  /** @brief Binds active layers in authored order and interns distinct sources.
   * @param current Current palette.
   * @param previous Retained palette.
   * @param t Current endpoint weight.
   * @param same Whether both palette endpoints have equal entries.
   * @param layers Consecutive canonical modulation stages.
   * @pre Every selected palette and active spatial source is prepared. */
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
        simple_palette(identity_controls &&
                       (!blend || (same && same_hue_table))),
        rgb_weight(blend && !same ? frac_to_q16(t) : 0) {
    const bool bound = bind_layers(layers);
    HS_CHECK(bound, "invalid ordered modulation column");
  }

  template <size_t N> void bind_specialization() {
    if constexpr (N == 3)
      hue_lift_mask = bind_hue_lift_mask();
  }

  template <size_t N>
  HS_HOT_INLINE Color4 shade_bound(const FieldSample &field) const {
    if constexpr (N == 3)
      if (hue_lift_mask)
        return shade_hue_lift_mask(field);
    return shade_ordered(field);
  }
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
