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

/// Scalar storage expected by a concrete source provider.
enum class LayerSampling : uint8_t {
  DYNAMIC,
  DIRECT,
  DIRECT_OR_PATH,
  BYTE_CUBE
};

/// Source storage and mechanic in a generated preset.
struct LayerLayout {
  /// Sampling strategy of the compiled layer.
  LayerSampling source = LayerSampling::DYNAMIC;
  /// Color operation of the compiled layer.
  ColorModulation::Mechanic mechanic = ColorModulation::Mechanic::HUE_ROTATION;
};

namespace StackDetail {
/** @brief Identifies the effective interpolation data of two borrowed cubes.
 * @param a First cube view.
 * @param b Second cube view.
 * @return Whether the cubes share the same interpolation data. */
inline bool same_cube(const Pattern::SphereLutView &a,
                      const Pattern::SphereLutView &b) {
  return a.linear_data == b.linear_data && (a.linear_data || a.data == b.data);
}

/** @brief Obtains the additional state of an exact phase representation.
 * @param phase Held phase representation.
 * @return Pointer to the expected phase type, or null. */
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
  enum class SourceKind : uint8_t {
    DIRECT,
    TABLE_TAP,
    BYTE_CUBE,
    POSITIVE_CUBE,
    FLOAT_CUBE
  };
  struct Source {
    union {
      const Modulation::Prepared *value;
      const int8_t *bytes;
      const float *linear;
    };
    SourceKind kind;
    uint8_t peer;
    bool positive_only;
  };
  struct Operation {
    ColorModulation::Params mechanic;
    uint8_t current;
    uint8_t previous;
    const PathModulationMix *path_mix;
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
  bool uses_path = false;
  bool varying_hue = false;
  bool varying_lift = false;
  bool varying_alpha = false;
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

  /** @brief Normalizes a sphere direction, including binary32 values close to unit length. */
  static HS_HOT_INLINE Pattern::SphereDirection
  sample_direction(const math::Vector &v) {
    const float squared = v.x * v.x + v.y * v.y + v.z * v.z;
    constexpr uint32_t ONE = std::bit_cast<uint32_t>(1.0f);
    constexpr uint32_t RANGE = 32;
    const uint32_t bits = std::bit_cast<uint32_t>(squared);
    const uint32_t slot = bits - (ONE - RANGE);
    if (slot >= 2 * RANGE)
      return Pattern::SphereDirection::admit(v);
    // In this range the negative quadratic term in sqrt(1+d) selects
    // the lower midpoint on odd bit offsets and leaves even offsets exact
    // after rounding.
    const uint32_t root_bits = (bits >> 1) + (ONE >> 1);
    if (root_bits == ONE)
      return Pattern::SphereDirection::from_unit(v);
    if (root_bits == ONE - 1) {
      // Dividing a normal binary32 by this root increments its magnitude
      // by one representable step. Subnormals round unchanged.
      const auto divide = [](float component) HS_HOT_INLINE {
        const auto bits = std::bit_cast<uint32_t>(component);
        return std::bit_cast<float>(bits + ((bits & 0x7f800000u) != 0));
      };
      return Pattern::SphereDirection::from_unit(
          {divide(v.x), divide(v.y), divide(v.z)});
    }
#if defined(__arm__) && defined(__RECIPROCAL_MATH__)
    static constexpr auto scale = [] {
      std::array<float, RANGE> result{};
      for (uint32_t i = 0; i < RANGE; ++i)
        result[i] = 1.0f / std::bit_cast<float>(ONE - RANGE / 2 + i);
      return result;
    }();
    const float inverse = scale[slot >> 1];
    return Pattern::SphereDirection::from_unit(
        {v.x * inverse, v.y * inverse, v.z * inverse});
#else
    const float root = std::bit_cast<float>(root_bits);
    return Pattern::SphereDirection::from_unit(
        {v.x / root, v.y / root, v.z / root});
#endif
  }

  static SourceKind kind(const Modulation::Prepared &endpoint) {
    if (endpoint.evaluate != Modulation::sample_cube || !endpoint.cube.active)
      return endpoint.immutable_table &&
                     endpoint.immutable_table->evaluator == endpoint.evaluate &&
                     endpoint.immutable_table->sample_tap
                 ? SourceKind::TABLE_TAP
                 : SourceKind::DIRECT;
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
      if ((type == SourceKind::DIRECT || type == SourceKind::TABLE_TAP) &&
          Modulation::same_direct_sample(value, other))
        return i;
    }
    const uint8_t i = source_count++;
    HS_AUDIT_CHECK(i < SOURCE_CAPACITY);
    sources[i] = {{&value}, type, peer, true};
    uses_cube |= type != SourceKind::DIRECT;
    return i;
  }

  /** @brief Resolves direct table pointers after every consumer has been interned. */
  void finish_sources() {
    for (uint8_t i = 0; i < source_count; ++i) {
      auto &source = sources[i];
      if (source.kind == SourceKind::BYTE_CUBE) {
        source.bytes = source.value->cube.data;
        if (source.positive_only)
          source.kind = SourceKind::POSITIVE_CUBE;
        if (source.peer != PATH) {
          const auto peer = sources[source.peer].kind;
          if (peer != SourceKind::BYTE_CUBE &&
              !(source.positive_only && peer == SourceKind::POSITIVE_CUBE))
            source.peer = PATH;
        }
      } else if (source.kind == SourceKind::FLOAT_CUBE)
        source.linear = source.value->cube.linear_data;
    }
  }

  static HS_HOT_INLINE float sample_bytes(const Pattern::SphereLutSample &tap,
                                          const int8_t *cube) {
    return Pattern::bounded_scalar(tap.sample_data(cube));
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

  HS_HOT_FLASH_MEMBER Color4 finish_simple(const FieldSample &field, float hue,
                                           float lift, float alpha) const {
    if (lift == 0.0f && hue != 0.0f && base.hue_rotation.active) {
      Pixel rotated;
      if (sample_guarded_hue_lut(base.hue_rotation.data,
                                 base.hue_rotation.gamut_mask, field.value, hue,
                                 rotated))
        return Color4(rotated, alpha);
    }
    const Pixel rgb = base.palette->get_color(field.value);
    if (lift == 0.0f && hue == 0.0f)
      return Color4(rgb, alpha);
    return Color4(finish_hue_lift_pixel(rgb, nullptr, lift, hue, 1.0f), alpha);
  }

  HS_HOT_FLASH_MEMBER Color4 shade_ordered(const FieldSample &field) const {
    float previous_alpha;
    ColorSample result;
    if (simple_palette) {
      previous_alpha =
          (constant_alpha >= 0.0f ? constant_alpha
                                  : base.palette->get_alpha(field.value)) *
          field.coverage;
      result.color.alpha = previous_alpha;
    } else
      result = palette_base(field, previous_alpha);
    std::array<float, SOURCE_CAPACITY + 1> values;
    values[PATH] = field.path_length;
    if (source_count) {
      const auto direction = sample_direction(field.sphere);
      Pattern::SphereLutSample tap{};
      if (uses_cube)
        tap = Pattern::prepare_sphere_lut_sample(direction);
      for (uint8_t i = 0; i < source_count; ++i) {
        const auto &source = sources[i];
        switch (source.kind) {
        case SourceKind::DIRECT:
          values[i] = Modulation::sample(*source.value, direction);
          break;
        case SourceKind::TABLE_TAP:
          values[i] = Pattern::bounded_scalar(
              source.value->immutable_table->sample_tap(*source.value, tap));
          break;
        case SourceKind::BYTE_CUBE:
        case SourceKind::POSITIVE_CUBE:
          if (source.peer != PATH &&
              tap.same_byte_corners(source.bytes, sources[source.peer].bytes))
            values[i] = values[source.peer];
          else
            values[i] = source.kind == SourceKind::POSITIVE_CUBE
                            ? positive_bytes(tap, source.bytes)
                            : sample_bytes(tap, source.bytes);
          break;
        case SourceKind::FLOAT_CUBE:
          values[i] = Pattern::bounded_scalar(tap.sample_data(source.linear));
          break;
        }
      }
    }
    float previous_hue = 0.0f;
    float previous_lift = 0.0f;
    for (uint8_t i = 0; i < op_count; ++i) {
      const auto &op = ops[i];
      const auto read = [&](uint8_t index) {
        return op.path_mix
                   ? op.path_mix->sample(op.mechanic.kind, values[index],
                                         field.path_length)
                   : values[index];
      };
      const float scalar = read(op.current);
      switch (op.mechanic.kind) {
      case ColorModulation::Mechanic::HUE_ROTATION:
        result.pending_hue =
            accumulate_hue(result.pending_hue, scalar, op.mechanic);
        if (varying_hue)
          previous_hue =
              accumulate_hue(previous_hue, read(op.previous), op.mechanic);
        break;
      case ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT:
        result.pending_lift =
            accumulate_lift(result.pending_lift, scalar, op.mechanic);
        if (varying_lift)
          previous_lift =
              accumulate_lift(previous_lift, read(op.previous), op.mechanic);
        break;
      case ColorModulation::Mechanic::OPACITY:
        result.color.alpha =
            opacity(result.color.alpha, scalar, op.mechanic.amount);
        if (varying_alpha)
          previous_alpha =
              opacity(previous_alpha, read(op.previous), op.mechanic.amount);
        break;
      }
    }
    if (varying_hue)
      result.pending_hue =
          blend_unit_hue(previous_hue, result.pending_hue, fraction);
    if (varying_lift)
      result.pending_lift =
          hs::lerp(previous_lift, result.pending_lift, fraction);
    if (varying_alpha)
      result.color.alpha =
          hs::lerp(previous_alpha, result.color.alpha, fraction);
    if (simple_palette)
      return finish_simple(field, result.pending_hue, result.pending_lift,
                           result.color.alpha);
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
      } else
        uses_path = true;
      using Mechanic = ColorModulation::Mechanic;
      const auto mechanic = layer->mechanic.kind;
      if (current != PATH && mechanic != Mechanic::POSITIVE_LIGHTNESS_LIFT) {
        sources[current].positive_only = false;
        sources[previous].positive_only = false;
      }
      if (previous != current) {
        varying_hue |= mechanic == Mechanic::HUE_ROTATION;
        varying_lift |= mechanic == Mechanic::POSITIVE_LIGHTNESS_LIFT;
        varying_alpha |= mechanic == Mechanic::OPACITY;
      }
      uses_path |= layer->path_mix.active();
      ops[op_count++] = {layer->mechanic, current, previous,
                         layer->path_mix.active() ? &layer->path_mix : nullptr};
    }
    return true;
  }

public:
  /// Whether an active layer samples displacement.
  /// @return Whether any layer consumes displacement path length.
  bool reads_path_length() const { return uses_path; }

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
    varying_alpha |= blend && !same_palette;
    if (simple_palette)
      constant_alpha = base.palette->constant_alpha().value_or(-1.0f);
    finish_sources();
  }

  /** @brief Tests whether each declared source can be sampled at its own layer.
   * @param layout Compiled source and mechanic declarations in layer order.
   * @return Whether independent per-layer sampling is compatible. */
  HS_SIZE_FLASH_MEMBER bool
  bind_independent_sources(std::span<const LayerLayout> layout) {
    if (!simple_palette || !source_count || layout.size() != op_count)
      return false;
    unsigned seen = 0;
    for (size_t i = 0; i < layout.size(); ++i) {
      const auto &op = ops[i];
      if (op.mechanic.kind != layout[i].mechanic)
        return false;
      if ((op.path_mix || op.current == PATH) &&
          layout[i].source != LayerSampling::DIRECT_OR_PATH)
        return false;
      if (op.current == PATH)
        continue;
      const auto compatible = [&](uint8_t source) {
        const auto k = sources[source].kind;
        return layout[i].source != LayerSampling::BYTE_CUBE
                   ? k == SourceKind::DIRECT || k == SourceKind::TABLE_TAP
                   : k == SourceKind::BYTE_CUBE ||
                         k == SourceKind::POSITIVE_CUBE;
      };
      const unsigned consumers = (1u << op.current) | (1u << op.previous);
      if ((seen & consumers) || !compatible(op.current) ||
          !compatible(op.previous))
        return false;
      seen |= consumers;
    }
    // Independent sources are interned in layer order. Expand pairs backwards.
    for (size_t i = op_count; i > 0; --i) {
      const auto &op = ops[i - 1];
      if (op.current == PATH) {
        sources[2 * i - 2] = {{nullptr}, SourceKind::DIRECT, PATH, false};
        continue;
      }
      const auto current = sources[op.current];
      sources[2 * i - 1] = sources[op.previous];
      sources[2 * i - 2] = current;
    }
    return true;
  }

  /// Accumulated color adjustments for a time endpoint.
  struct Adjustments {
    /** @var Adjustments::lift
     * @brief Accumulated lightness lift. */
    /** @var Adjustments::alpha
     * @brief Accumulated opacity. */
    /// Hue rotation in turns.
    float hue = 0.0f, lift = 0.0f, alpha;
  };

  /** @brief Accumulates the selected layer at the two time endpoints.
   * @param op Bound modulation operation.
   * @param direction Normalized source direction.
   * @param tap Shared cube interpolation tap.
   * @param path_length Accumulated displacement path length.
   * @param current Incoming endpoint adjustments.
   * @param previous Outgoing endpoint adjustments. */
  template <LayerLayout Layout, bool Unique, size_t I>
  HS_HOT_INLINE void
  apply_source(const Operation &op, const Pattern::SphereDirection &direction,
               const Pattern::SphereLutSample &tap, float path_length,
               Adjustments &current, Adjustments &previous) const {
    const auto read = [&](uint8_t index) HS_HOT_INLINE {
      const auto &source = sources[index];
      if constexpr (Layout.source != LayerSampling::BYTE_CUBE) {
        if constexpr (Layout.source == LayerSampling::DIRECT_OR_PATH)
          if (!source.value)
            return path_length;
        const auto &prepared = *source.value;
        const float spatial =
            source.kind == SourceKind::TABLE_TAP
                ? Pattern::bounded_scalar(
                      prepared.immutable_table->sample_tap(prepared, tap))
                : Modulation::sample(prepared, direction);
        if constexpr (Layout.source == LayerSampling::DIRECT_OR_PATH)
          if (op.path_mix)
            return op.path_mix->sample<Layout.mechanic>(spatial, path_length);
        return spatial;
      } else {
        if constexpr (Layout.mechanic ==
                      ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT)
          return positive_bytes(tap, source.bytes);
        else
          return sample_bytes(tap, source.bytes);
      }
    };
    using M = ColorModulation::Mechanic;
    const auto accumulate = [&](float prior, float scalar) HS_HOT_INLINE {
      if constexpr (Layout.mechanic == M::HUE_ROTATION) {
        if constexpr (Unique)
          return math::wrap_t(op.mechanic.amount * scalar);
        return accumulate_hue(prior, scalar, op.mechanic);
      } else if constexpr (Layout.mechanic == M::POSITIVE_LIGHTNESS_LIFT) {
        const float positive = Layout.source == LayerSampling::BYTE_CUBE
                                   ? scalar
                                   : hs::clamp(scalar, 0.0f, 1.0f);
        const float lift = op.mechanic.amount * positive;
        if constexpr (Unique)
          return lift;
        return prior + (1.0f - prior) * lift;
      } else if constexpr (Layout.source == LayerSampling::DIRECT_OR_PATH)
        return opacity(prior, scalar, op.mechanic.amount);
      else
        return prior *
               hs::lerp(1.0f, 0.5f * (scalar + 1.0f), op.mechanic.amount);
    };
    constexpr auto member =
        Layout.mechanic == M::HUE_ROTATION ? &Adjustments::hue
        : Layout.mechanic == M::OPACITY    ? &Adjustments::alpha
                                           : &Adjustments::lift;
    const float scalar = read(2 * I);
    const float value = accumulate(current.*member, scalar);
    if constexpr (Unique) {
      bool same;
      if constexpr (Layout.source == LayerSampling::BYTE_CUBE) {
        const auto *current_bytes = sources[2 * I].bytes;
        const auto *previous_bytes = sources[2 * I + 1].bytes;
        same = current_bytes == previous_bytes ||
               tap.same_byte_corners(current_bytes, previous_bytes);
      } else {
        same = op.previous == op.current;
      }
      if (same)
        current.*member = value;
      else {
        const float old = accumulate(current.*member, read(2 * I + 1));
        if constexpr (Layout.mechanic == M::HUE_ROTATION)
          current.*member = blend_unit_hue(old, value, fraction);
        else
          current.*member = hs::lerp(old, value, fraction);
      }
    } else {
      current.*member = value;
      const float old = op.previous == op.current ? scalar : read(2 * I + 1);
      previous.*member = accumulate(previous.*member, old);
    }
  }

  /// Counts the compiled layers that use the selected mechanic.
  template <ColorModulation::Mechanic Mechanic, LayerLayout... Layout>
  static constexpr size_t uses =
      (size_t{0} + ... + size_t(Layout.mechanic == Mechanic));

  /** @brief Samples each independent source using its provider's storage layout.
   * @param field Scalar field and sample provenance.
   * @return Color from independently sampled layers. */
  template <LayerLayout... Layout, size_t... I>
  HS_HOT_ITCM_MEMBER Color4 shade_independent(const FieldSample &field,
                                              std::index_sequence<I...>) const {
    const auto direction = sample_direction(field.sphere);
    Pattern::SphereLutSample tap{};
    if constexpr (((Layout.source == LayerSampling::BYTE_CUBE) || ...))
      tap = Pattern::prepare_sphere_lut_sample(direction);
    else if (uses_cube)
      tap = Pattern::prepare_sphere_lut_sample(direction);
    const float alpha =
        (constant_alpha >= 0.0f ? constant_alpha
                                : base.palette->get_alpha(field.value)) *
        field.coverage;
    Adjustments current{0.0f, 0.0f, alpha}, previous = current;
    (apply_source<Layout, uses<Layout.mechanic, Layout...> == 1, I>(
         ops[I], direction, tap, field.path_length, current, previous),
     ...);
    if constexpr (uses<ColorModulation::Mechanic::HUE_ROTATION, Layout...> > 1)
      if (varying_hue)
        current.hue = blend_unit_hue(previous.hue, current.hue, fraction);
    if constexpr (uses<ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT,
                       Layout...> > 1)
      if (varying_lift)
        current.lift = hs::lerp(previous.lift, current.lift, fraction);
    if constexpr (uses<ColorModulation::Mechanic::OPACITY, Layout...> > 1)
      if (varying_alpha)
        current.alpha = hs::lerp(previous.alpha, current.alpha, fraction);
    return finish_simple(field, current.hue, current.lift, current.alpha);
  }

  /** @brief Shades with the ordered layers bound for this frame.
   * @param field Scalar field and sample provenance.
   * @return Color after modulation and brightness. */
  HS_HOT_INLINE Color4 shade_bound(const FieldSample &field) const {
    return shade_ordered(field);
  }
};

/** @brief Column binding for a typed sequence of modulation stages.
 * @tparam N Number of consecutive stages. */
template <size_t N, LayerLayout... Layout>
struct ModulationStack : OrderedModulationKernel {
private:
  static constexpr bool typed =
      sizeof...(Layout) == N &&
      ((Layout.source != LayerSampling::DYNAMIC) && ...);
  bool independent = [&] {
    if constexpr (typed) {
      static constexpr std::array layout{Layout...};
      return bind_independent_sources(layout);
    } else
      return false;
  }();

public:
  /** @brief Binds the ordered modulation stages to a generated palette.
   * @param palette Prepared temporal palette.
   * @param same_palette Whether the palette entries match.
   * @param layers Prepared layer endpoints. */
  ModulationStack(const TemporalGeneratedPaletteState &palette,
                  bool same_palette,
                  std::array<const ModulationState *, N> layers)
      : OrderedModulationKernel(palette.current, palette.previous,
                                palette.fraction, same_palette, layers) {}

  /** @brief Binds consecutive layers on one non-temporal palette.
   * @param palette Current palette.
   * @param layers Prepared layer endpoints. */
  ModulationStack(const GeneratedPaletteState &palette,
                  std::array<const ModulationState *, N> layers)
      : OrderedModulationKernel(palette, palette, 1.0f, true, layers) {}

  /** @brief Evaluates the bound column at one field sample.
   * @param field Shared source direction, coverage and displacement.
   * @return Finished straight-alpha color. */
  HS_HOT_INLINE Color4 shade(const FieldSample &field) const {
    if constexpr (typed)
      if (independent)
        return shade_independent<Layout...>(field,
                                            std::make_index_sequence<N>{});
    return shade_bound(field);
  }
};

} // namespace Pullback::Color
