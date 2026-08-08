/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "color/color.h"

/**
 * @brief A compiled perceptual palette.
 */
class GenerativePalette : public Palette {
public:
  struct ControlKey {
    float L;
    union {
      float chroma;
      float C;
      float q;
    };
    float h;

    constexpr ControlKey(float lightness = 0.0f, float chroma = 0.0f,
                         float hue = 0.0f)
        : L(lightness), chroma(chroma), h(hue) {}
  };

  struct Snapshot {
    /** Compact owned state that keeps palette animations inside their budget. */
    struct Key {
      std::array<uint8_t, 7> bytes{};
    };

    std::array<Key, PALETTE_MAX_KEYS> keys;
    float lightness_low;
    float lightness_high;
    float chroma_low;
    float chroma_high;
    uint8_t key_count;
    AxisCurve lightness_curve;
    AxisCurve chroma_curve;
    uint8_t reserved = 0;
  };

  struct Diagnostic {
    float L;
    float C;
    float q;
    float C_max;
    float h_path;
    float h_final;
    bool fallback_mapped;
  };

  static_assert(sizeof(Snapshot) == 48);

  HS_COLD_MEMBER GenerativePalette() { initialize(PaletteRecipe{}); }

  HS_COLD_MEMBER explicit GenerativePalette(const PaletteRecipe &recipe) {
    PaletteRecipe canonical;
    PaletteCompileStatus status;
    const bool compiled = canonicalize(recipe, canonical, status);
    HS_CHECK(compiled, "GenerativePalette recipe error %d at field %d",
             static_cast<int>(status.code), static_cast<int>(status.field));
    initialize(canonical);
  }

  static bool try_compile(const PaletteRecipe &input, GenerativePalette &output,
                          PaletteRecipe &canonical,
                          PaletteCompileStatus &status) {
    PaletteRecipe next_recipe;
    if (!canonicalize(input, next_recipe, status))
      return false;

    GenerativePalette next(Unchecked{}, next_recipe);
    output = next;
    canonical = next_recipe;
    return true;
  }

  Snapshot snapshot() const {
    Snapshot result{};
    result.key_count = key_count;
    result.lightness_low = lightness_axis.low;
    result.lightness_high = lightness_axis.high;
    result.chroma_low = chroma_axis.low;
    result.chroma_high = chroma_axis.high;
    result.lightness_curve = lightness_axis.curve;
    result.chroma_curve = chroma_axis.curve;
    for (int i = 0; i < key_count; ++i)
      result.keys[i] = encode_snapshot_key(keys[i]);
    return result;
  }

  static ControlKey snapshot_key(const Snapshot &snapshot, int index) {
    index = hs::clamp(index, 0, static_cast<int>(snapshot.key_count) - 1);
    ControlKey key = decode_snapshot_key(snapshot.keys[index]);
    const float position = index / float(snapshot.key_count - 1);
    if (snapshot.lightness_curve != AxisCurve::CUSTOM)
      key.L = evaluate_axis(snapshot.lightness_low, snapshot.lightness_high,
                            snapshot.lightness_curve, position);
    if (snapshot.chroma_curve != AxisCurve::CUSTOM)
      key.chroma = evaluate_axis(snapshot.chroma_low, snapshot.chroma_high,
                                 snapshot.chroma_curve, position);
    return key;
  }

  PaletteDomain palette_domain() const { return domain; }
  ChromaBasis palette_chroma_basis() const { return chroma_basis; }
  ColorPath palette_color_path() const { return color_path; }
  SegmentEase segment_easing() const { return easing; }
  float palette_headroom() const { return headroom; }
  float palette_hue_torsion() const { return hue_torsion; }
  float palette_input_offset() const { return input_offset; }
  float palette_input_span() const { return input_span; }
  uint8_t palette_key_count() const { return key_count; }

  bool mirrors_domain() const { return domain == PaletteDomain::MIRROR; }
  bool loops_domain() const { return domain == PaletteDomain::LOOP; }

  Diagnostic diagnose(float t) const {
    const Evaluated value = finish_evaluation(evaluate_path(t));
    return {value.lab.L,          value.C,      value.q,
            value.C_max,          value.h_path, value.h_final,
            value.fallback_mapped};
  }

  OKLCH resolved_oklch_key(int index) const {
    index = hs::clamp(index, 0, static_cast<int>(key_count) - 1);
    ControlKey key = keys[index];
    const float position = index / float(key_count - 1);
    if (lightness_axis.curve != AxisCurve::CUSTOM)
      key.L = evaluate_axis(lightness_axis.low, lightness_axis.high,
                            lightness_axis.curve, position);
    if (chroma_axis.curve != AxisCurve::CUSTOM)
      key.chroma = evaluate_axis(chroma_axis.low, chroma_axis.high,
                                 chroma_axis.curve, position);
    const float h_final = key.h + hue_torsion * (key.L - 0.5f);
    return {key.L, realized_chroma(key, h_final), h_final};
  }

  /**
   * @brief True when this palette and @p other can be key-morphed by lerp().
   * @details Requires identical evaluation policy — domain, easing, color
   * path, chroma basis, complementary evaluation, key count, headroom,
   * torsion, falloff, input window — plus corresponding segment hue deltas
   * within half a turn of each other, so an interpolated arc never crosses a
   * whole-turn ambiguity or collapses through zero. A loop additionally
   * requires the same integer closing travel, or its seam breaks mid-morph.
   * Incompatible palettes must transition through a baked crossfade
   * (BakedPalette::rebake_crossfade) instead.
   */
  bool morph_compatible(const GenerativePalette &other) const {
    if (domain != other.domain || easing != other.easing ||
        color_path != other.color_path || chroma_basis != other.chroma_basis ||
        complementary_harmony != other.complementary_harmony ||
        key_count != other.key_count || headroom != other.headroom ||
        hue_torsion != other.hue_torsion ||
        falloff_start != other.falloff_start ||
        input_offset != other.input_offset || input_span != other.input_span)
      return false;
    for (int i = 1; i < key_count; ++i) {
      if (fabsf((keys[i].h - keys[i - 1].h) -
                (other.keys[i].h - other.keys[i - 1].h)) >= PI_F)
        return false;
    }
    if (domain == PaletteDomain::LOOP) {
      if (fabsf((closing_hue - keys[key_count - 1].h) -
                (other.closing_hue - other.keys[key_count - 1].h)) >= PI_F)
        return false;
      if (loop_turns() != other.loop_turns())
        return false;
    }
    return true;
  }

  /**
   * @brief Morphs this palette between two morph-compatible palettes.
   * @param from Source palette (amount = 0).
   * @param to Target palette (amount = 1).
   * @param amount Blend amount; clamped to [0, 1].
   * @details Adopts @p from's evaluation policy, snapshot-lerps the keys, and
   * interpolates the loop-closing hue delta, which snapshots do not carry.
   * Callers gate on morph_compatible().
   */
  HS_COLD_MEMBER void lerp(const GenerativePalette &from,
                           const GenerativePalette &to, float amount) {
    amount = hs::clamp(amount, 0.0f, 1.0f);
    if (amount == 0.0f) {
      *this = from;
      return;
    }
    if (amount == 1.0f) {
      *this = to;
      return;
    }
    const float closing_from = from.closing_hue - from.keys[0].h;
    const float closing_to = to.closing_hue - to.keys[0].h;
    *this = from;
    lerp(from.snapshot(), to.snapshot(), amount);
    closing_hue =
        keys[0].h + closing_from + (closing_to - closing_from) * amount;
  }

  HS_COLD_MEMBER void lerp(const Snapshot &from, const Snapshot &to,
                           float amount) {
    amount = hs::clamp(amount, 0.0f, 1.0f);
    if (amount == 0.0f) {
      assign(from);
      return;
    }
    if (amount == 1.0f) {
      assign(to);
      return;
    }

    // Mixed axis curves morph through CUSTOM: the per-key samples below carry
    // each side's own curve, while a non-CUSTOM curve would ignore them and
    // evaluate the lerped low/high — wrong for the CUSTOM side, and a pop at
    // the endpoint copy.
    lightness_axis = {
        from.lightness_low + (to.lightness_low - from.lightness_low) * amount,
        from.lightness_high +
            (to.lightness_high - from.lightness_high) * amount,
        from.lightness_curve == to.lightness_curve ? from.lightness_curve
                                                   : AxisCurve::CUSTOM,
    };
    chroma_axis = {
        from.chroma_low + (to.chroma_low - from.chroma_low) * amount,
        from.chroma_high + (to.chroma_high - from.chroma_high) * amount,
        from.chroma_curve == to.chroma_curve ? from.chroma_curve
                                             : AxisCurve::CUSTOM,
    };
    key_count = std::max(from.key_count, to.key_count);
    std::array<ControlKey, PALETTE_MAX_KEYS> from_keys{};
    std::array<ControlKey, PALETTE_MAX_KEYS> to_keys{};
    for (int i = 0; i < key_count; ++i) {
      const float position = key_count > 1 ? i / float(key_count - 1) : 0.0f;
      from_keys[i] = sample_snapshot(from, position);
      to_keys[i] = sample_snapshot(to, position);
    }
    int anchor = -1;
    for (int i = 0; i < key_count; ++i) {
      if (is_chromatic(from_keys[i]) && is_chromatic(to_keys[i])) {
        anchor = i;
        break;
      }
    }

    const float anchor_delta =
        anchor >= 0 ? wrap_angle_pi(to_keys[anchor].h - from_keys[anchor].h)
                    : 0.0f;
    for (int i = 0; i < key_count; ++i) {
      const bool from_chromatic = is_chromatic(from_keys[i]);
      const bool to_chromatic = is_chromatic(to_keys[i]);
      float hue = 0.0f;
      if (from_chromatic && to_chromatic) {
        const float delta =
            anchor < 0
                ? wrap_angle_pi(to_keys[i].h - from_keys[i].h)
                : anchor_delta +
                      wrap_angle_pi((to_keys[i].h - to_keys[anchor].h) -
                                    (from_keys[i].h - from_keys[anchor].h));
        hue = from_keys[i].h + delta * amount;
      } else if (from_chromatic) {
        hue = from_keys[i].h;
      } else if (to_chromatic) {
        hue = to_keys[i].h;
      }
      keys[i] = {from_keys[i].L + (to_keys[i].L - from_keys[i].L) * amount,
                 from_keys[i].chroma +
                     (to_keys[i].chroma - from_keys[i].chroma) * amount,
                 hue};
    }
    for (int i = key_count; i < PALETTE_MAX_KEYS; ++i)
      keys[i] = {};
  }

  HS_COLD_MEMBER Color4 get(float t) const override {
    const LinRGB rgb = oklab_to_linear_rgb_gamut(evaluate_path(t).lab);
    return Color4(Pixel(float_to_pixel16(rgb.r), float_to_pixel16(rgb.g),
                        float_to_pixel16(rgb.b)),
                  1.0f);
  }

private:
  struct Unchecked {};

  struct Segment {
    ControlKey left;
    ControlKey right;
    ControlKey envelope;
    float visibility;
    float progress;
    float relationship_position;
  };

  struct Evaluated {
    OKLab lab;
    float C;
    float q;
    float C_max;
    float h_path;
    float h_final;
    bool fallback_mapped;
  };

  struct PathEvaluation {
    OKLab lab;
    float q;
    float h_path;
    float h_final;
    float C_max;
  };

  struct AxisState {
    float low;
    float high;
    AxisCurve curve;
  };

  static constexpr float MAX_SWEEP_TURNS = 16.0f;
  static constexpr float MAX_CUSTOM_DELTA = 16.0f;
  static constexpr float MAX_CUSTOM_ABS_INPUT = 4096.0f;
  static constexpr float MAX_ABSOLUTE_CHROMA = 1.0f;
  static constexpr float MAX_ABS_TORSION = 4.0f * PI_F;
  static constexpr float SWEEP_INTEGER_EPS = 1e-6f;
  static constexpr float SNAPSHOT_AXIS_STEPS = 4095.0f;

  static float evaluate_axis(float low, float high, AxisCurve curve,
                             float position) {
    position = hs::clamp(position, 0.0f, 1.0f);
    switch (curve) {
    case AxisCurve::CONSTANT:
      return low;
    case AxisCurve::ASCENDING:
      return low + (high - low) * position;
    case AxisCurve::DESCENDING:
      return high + (low - high) * position;
    case AxisCurve::BELL:
      return low + (high - low) * (1.0f - fabsf(position * 2.0f - 1.0f));
    case AxisCurve::CUP:
      return high - (high - low) * (1.0f - fabsf(position * 2.0f - 1.0f));
    case AxisCurve::CUSTOM:
      return 0.0f;
    }
    return 0.0f;
  }

  static AxisState axis_state(const AxisControls &controls) {
    if (controls.curve == AxisCurve::CONSTANT) {
      const float value = hs::clamp(controls.center, 0.0f, 1.0f);
      return {value, value, controls.curve};
    }
    return {hs::clamp(controls.center - controls.range * 0.5f, 0.0f, 1.0f),
            hs::clamp(controls.center + controls.range * 0.5f, 0.0f, 1.0f),
            controls.curve};
  }

  static AxisState axis_state(const ChromaControls &controls) {
    AxisControls axis;
    axis.curve = controls.curve;
    axis.center = controls.center;
    axis.range = controls.range;
    return axis_state(axis);
  }

  HS_COLD_MEMBER GenerativePalette(Unchecked, const PaletteRecipe &recipe) {
    initialize(recipe);
  }

  static uint64_t field_bit(PaletteRecipeField field) {
    return uint64_t{1} << static_cast<uint8_t>(field);
  }

  static bool fail(PaletteCompileStatus &status, PaletteCompileCode code,
                   PaletteRecipeField field) {
    status.code = code;
    status.field = field;
    return false;
  }

  template <typename Enum> static bool enum_in_range(Enum value, Enum last) {
    return static_cast<uint8_t>(value) <= static_cast<uint8_t>(last);
  }

  static bool finite(float value, PaletteRecipeField field,
                     PaletteCompileStatus &status) {
    return std::isfinite(value) ||
           fail(status, PaletteCompileCode::NON_FINITE, field);
  }

  HS_COLD_MEMBER static void clamp_field(float &value, float low, float high,
                                         PaletteRecipeField field,
                                         PaletteCompileStatus &status) {
    const float clamped = hs::clamp(value, low, high);
    if (clamped != value) {
      value = clamped;
      status.adjustments.clamped_fields |= field_bit(field);
    }
  }

  static bool validate_enums(const PaletteRecipe &recipe,
                             PaletteCompileStatus &status) {
    if (!enum_in_range(recipe.domain, PaletteDomain::LOOP))
      return fail(status, PaletteCompileCode::INVALID_ENUM,
                  PaletteRecipeField::PALETTE_DOMAIN);
    if (!enum_in_range(recipe.easing, SegmentEase::SMOOTHSTEP))
      return fail(status, PaletteCompileCode::INVALID_ENUM,
                  PaletteRecipeField::EASING);
    if (!enum_in_range(recipe.color_path, ColorPath::OKLAB_CARTESIAN))
      return fail(status, PaletteCompileCode::INVALID_ENUM,
                  PaletteRecipeField::COLOR_PATH);
    if (!enum_in_range(recipe.hue.mode, HueMode::CUSTOM))
      return fail(status, PaletteCompileCode::INVALID_ENUM,
                  PaletteRecipeField::HUE_MODE);
    if (!enum_in_range(recipe.hue.harmony, PaletteHarmony::SQUARE))
      return fail(status, PaletteCompileCode::INVALID_ENUM,
                  PaletteRecipeField::HARMONY);
    if (!enum_in_range(recipe.hue.direction, HueDirection::COUNTERCLOCKWISE))
      return fail(status, PaletteCompileCode::INVALID_ENUM,
                  PaletteRecipeField::HUE_DIRECTION);
    if (!enum_in_range(recipe.lightness.curve, AxisCurve::CUSTOM))
      return fail(status, PaletteCompileCode::INVALID_ENUM,
                  PaletteRecipeField::LIGHTNESS_CURVE);
    if (!enum_in_range(recipe.chroma.curve, AxisCurve::CUSTOM))
      return fail(status, PaletteCompileCode::INVALID_ENUM,
                  PaletteRecipeField::CHROMA_CURVE);
    if (!enum_in_range(recipe.chroma.basis, ChromaBasis::ABSOLUTE))
      return fail(status, PaletteCompileCode::INVALID_ENUM,
                  PaletteRecipeField::CHROMA_BASIS);
    return true;
  }

  static uint8_t control_key_count(const PaletteRecipe &recipe) {
    if (recipe.hue.mode == HueMode::SWEEP)
      return 2;
    if (recipe.hue.mode == HueMode::CUSTOM)
      return 3;
    switch (recipe.hue.harmony) {
    case PaletteHarmony::ANALOGOUS:
    case PaletteHarmony::SPLIT_COMPLEMENTARY:
    case PaletteHarmony::TRIADIC:
      return 3;
    case PaletteHarmony::ACCENTED_ANALOGOUS:
    case PaletteHarmony::TETRADIC:
    case PaletteHarmony::SQUARE:
      return 4;
    case PaletteHarmony::MONOCHROMATIC:
    case PaletteHarmony::COMPLEMENTARY:
      return 2;
    }
    return 2;
  }

  static bool validate_finite(const PaletteRecipe &recipe,
                              PaletteCompileStatus &status) {
    if (!finite(recipe.hue.base_turns, PaletteRecipeField::BASE_TURNS,
                status) ||
        !finite(recipe.hue.spread_turns, PaletteRecipeField::SPREAD_TURNS,
                status) ||
        !finite(recipe.hue.sweep_turns, PaletteRecipeField::SWEEP_TURNS,
                status))
      return false;
    for (int i = 0; i < PALETTE_MAX_KEYS; ++i) {
      if (!finite(
              recipe.hue.custom_turns[i],
              static_cast<PaletteRecipeField>(
                  static_cast<uint8_t>(PaletteRecipeField::CUSTOM_TURNS_0) + i),
              status))
        return false;
    }
    if (!finite(recipe.lightness.center, PaletteRecipeField::LIGHTNESS_CENTER,
                status) ||
        !finite(recipe.lightness.range, PaletteRecipeField::LIGHTNESS_RANGE,
                status))
      return false;
    for (int i = 0; i < PALETTE_MAX_KEYS; ++i) {
      if (!finite(
              recipe.lightness.custom[i],
              static_cast<PaletteRecipeField>(
                  static_cast<uint8_t>(PaletteRecipeField::LIGHTNESS_CUSTOM_0) +
                  i),
              status))
        return false;
    }
    if (!finite(recipe.chroma.center, PaletteRecipeField::CHROMA_CENTER,
                status) ||
        !finite(recipe.chroma.range, PaletteRecipeField::CHROMA_RANGE, status))
      return false;
    for (int i = 0; i < PALETTE_MAX_KEYS; ++i) {
      if (!finite(
              recipe.chroma.custom[i],
              static_cast<PaletteRecipeField>(
                  static_cast<uint8_t>(PaletteRecipeField::CHROMA_CUSTOM_0) +
                  i),
              status))
        return false;
    }
    return finite(recipe.chroma.headroom, PaletteRecipeField::CHROMA_HEADROOM,
                  status) &&
           finite(recipe.hue_torsion, PaletteRecipeField::HUE_TORSION,
                  status) &&
           finite(recipe.falloff_start, PaletteRecipeField::FALLOFF_START,
                  status) &&
           finite(recipe.input.offset, PaletteRecipeField::INPUT_OFFSET,
                  status) &&
           finite(recipe.input.span, PaletteRecipeField::INPUT_SPAN, status);
  }

  HS_COLD_MEMBER static bool canonicalize(const PaletteRecipe &input,
                                          PaletteRecipe &canonical,
                                          PaletteCompileStatus &status) {
    status = {};
    if (input.schema_version != PaletteRecipe::SCHEMA_VERSION)
      return fail(status, PaletteCompileCode::INVALID_SCHEMA,
                  PaletteRecipeField::SCHEMA_VERSION);
    if (!validate_enums(input, status) || !validate_finite(input, status))
      return false;

    PaletteRecipe recipe = input;
    const uint8_t key_count = control_key_count(recipe);
    clamp_field(recipe.input.offset, 0.0f, 1.0f,
                PaletteRecipeField::INPUT_OFFSET, status);
    clamp_field(recipe.input.span, 0.0f, 1.0f - recipe.input.offset,
                PaletteRecipeField::INPUT_SPAN, status);
    const float wrapped_base =
        recipe.hue.base_turns - floorf(recipe.hue.base_turns);
    if (wrapped_base != recipe.hue.base_turns) {
      recipe.hue.base_turns = wrapped_base;
      status.adjustments.wrapped_fields |=
          field_bit(PaletteRecipeField::BASE_TURNS);
    }
    clamp_field(recipe.hue.spread_turns, 0.0f, 0.25f,
                PaletteRecipeField::SPREAD_TURNS, status);
    if (fabsf(recipe.hue.sweep_turns) > MAX_SWEEP_TURNS)
      return fail(status, PaletteCompileCode::HUE_LIMIT,
                  PaletteRecipeField::SWEEP_TURNS);

    for (int i = 0; i < PALETTE_MAX_KEYS; ++i) {
      if (fabsf(recipe.hue.custom_turns[i]) > MAX_CUSTOM_ABS_INPUT)
        return fail(
            status, PaletteCompileCode::HUE_LIMIT,
            static_cast<PaletteRecipeField>(
                static_cast<uint8_t>(PaletteRecipeField::CUSTOM_TURNS_0) + i));
    }
    const float custom_offset = floorf(recipe.hue.custom_turns[0]);
    if (custom_offset != 0.0f) {
      for (int i = 0; i < PALETTE_MAX_KEYS; ++i) {
        recipe.hue.custom_turns[i] -= custom_offset;
        status.adjustments.canonicalized_fields |=
            field_bit(static_cast<PaletteRecipeField>(
                static_cast<uint8_t>(PaletteRecipeField::CUSTOM_TURNS_0) + i));
      }
    }
    for (int i = 1; i < key_count; ++i) {
      if (fabsf(recipe.hue.custom_turns[i] - recipe.hue.custom_turns[i - 1]) >
          MAX_CUSTOM_DELTA)
        return fail(
            status, PaletteCompileCode::HUE_LIMIT,
            static_cast<PaletteRecipeField>(
                static_cast<uint8_t>(PaletteRecipeField::CUSTOM_TURNS_0) + i));
    }

    clamp_field(recipe.lightness.center, 0.0f, 1.0f,
                PaletteRecipeField::LIGHTNESS_CENTER, status);
    clamp_field(recipe.lightness.range, 0.0f, 1.0f,
                PaletteRecipeField::LIGHTNESS_RANGE, status);
    for (int i = 0; i < PALETTE_MAX_KEYS; ++i) {
      clamp_field(
          recipe.lightness.custom[i], 0.0f, 1.0f,
          static_cast<PaletteRecipeField>(
              static_cast<uint8_t>(PaletteRecipeField::LIGHTNESS_CUSTOM_0) + i),
          status);
    }

    const float chroma_max = recipe.chroma.basis == ChromaBasis::ABSOLUTE
                                 ? MAX_ABSOLUTE_CHROMA
                                 : 1.0f;
    clamp_field(recipe.chroma.center, 0.0f, chroma_max,
                PaletteRecipeField::CHROMA_CENTER, status);
    clamp_field(recipe.chroma.range, 0.0f, chroma_max,
                PaletteRecipeField::CHROMA_RANGE, status);
    for (int i = 0; i < PALETTE_MAX_KEYS; ++i) {
      clamp_field(
          recipe.chroma.custom[i], 0.0f, chroma_max,
          static_cast<PaletteRecipeField>(
              static_cast<uint8_t>(PaletteRecipeField::CHROMA_CUSTOM_0) + i),
          status);
    }
    if (recipe.chroma.basis == ChromaBasis::ABSOLUTE) {
      if (recipe.chroma.headroom != 1.0f) {
        recipe.chroma.headroom = 1.0f;
        status.adjustments.canonicalized_fields |=
            field_bit(PaletteRecipeField::CHROMA_HEADROOM);
      }
    } else {
      clamp_field(recipe.chroma.headroom, 0.0f, 1.0f,
                  PaletteRecipeField::CHROMA_HEADROOM, status);
    }

    if (fabsf(recipe.hue_torsion) > MAX_ABS_TORSION)
      return fail(status, PaletteCompileCode::HUE_LIMIT,
                  PaletteRecipeField::HUE_TORSION);
    if (recipe.domain == PaletteDomain::FALLOFF) {
      if (!(recipe.falloff_start > 2.0f / 3.0f && recipe.falloff_start < 1.0f))
        return fail(status, PaletteCompileCode::INVALID_FALLOFF_START,
                    PaletteRecipeField::FALLOFF_START);
    } else if (recipe.falloff_start != 0.90f) {
      recipe.falloff_start = 0.90f;
      status.adjustments.canonicalized_fields |=
          field_bit(PaletteRecipeField::FALLOFF_START);
    }

    if (recipe.chroma.basis == ChromaBasis::PATH_MINIMUM)
      return fail(status, PaletteCompileCode::INCOMPATIBLE_OPTIONS,
                  PaletteRecipeField::CHROMA_BASIS);
    if (recipe.domain == PaletteDomain::LOOP &&
        recipe.hue.mode == HueMode::SWEEP) {
      const float turns = roundf(recipe.hue.sweep_turns);
      if (fabsf(recipe.hue.sweep_turns - turns) > SWEEP_INTEGER_EPS)
        return fail(status, PaletteCompileCode::NON_INTEGER_LOOP_SWEEP,
                    PaletteRecipeField::SWEEP_TURNS);
      recipe.hue.sweep_turns = turns;
    }

    canonical = recipe;
    return true;
  }

  static float wrap01(float value) { return value - floorf(value); }

  /** @brief Whole turns of a loop's closing travel from the first key. */
  int loop_turns() const {
    return static_cast<int>(
        roundf((closing_hue - keys[0].h) * (1.0f / (2.0f * PI_F))));
  }

  HS_COLD_MEMBER static float directed_delta(float delta,
                                             HueDirection direction) {
    const float wrapped = wrap01(delta);
    switch (direction) {
    case HueDirection::SHORTEST:
      if (wrapped < 0.5f)
        return wrapped;
      if (wrapped > 0.5f)
        return wrapped - 1.0f;
      return delta < 0.0f ? -0.5f : 0.5f;
    case HueDirection::CLOCKWISE:
      return wrapped == 0.0f ? 0.0f : wrapped - 1.0f;
    case HueDirection::COUNTERCLOCKWISE:
      return wrapped;
    }
    return 0.0f;
  }

  HS_COLD_MEMBER static void resolve_axis(const AxisControls &controls,
                                          uint8_t count,
                                          float out[PALETTE_MAX_KEYS]) {
    const float low =
        hs::clamp(controls.center - controls.range * 0.5f, 0.0f, 1.0f);
    const float high =
        hs::clamp(controls.center + controls.range * 0.5f, 0.0f, 1.0f);
    for (int i = 0; i < count; ++i) {
      const float position = i / float(count - 1);
      switch (controls.curve) {
      case AxisCurve::CONSTANT:
        out[i] = controls.center;
        break;
      case AxisCurve::ASCENDING:
        out[i] = low + (high - low) * position;
        break;
      case AxisCurve::DESCENDING:
        out[i] = high + (low - high) * position;
        break;
      case AxisCurve::BELL:
        out[i] = low + (high - low) * (1.0f - fabsf(position * 2.0f - 1.0f));
        break;
      case AxisCurve::CUP:
        out[i] = high - (high - low) * (1.0f - fabsf(position * 2.0f - 1.0f));
        break;
      case AxisCurve::CUSTOM:
        out[i] = controls.custom[i];
        break;
      }
    }
  }

  HS_COLD_MEMBER static void resolve_axis(const ChromaControls &controls,
                                          uint8_t count,
                                          float out[PALETTE_MAX_KEYS]) {
    AxisControls axis;
    axis.curve = controls.curve;
    axis.center = controls.center;
    axis.range = controls.range;
    for (int i = 0; i < PALETTE_MAX_KEYS; ++i)
      axis.custom[i] = controls.custom[i];
    resolve_axis(axis, count, out);
  }

  HS_COLD_MEMBER static void resolve_harmony(const PaletteRecipe &recipe,
                                             float hues[PALETTE_MAX_KEYS]) {
    const float base = recipe.hue.base_turns;
    const float spread = recipe.hue.spread_turns;
    const float orientation =
        recipe.hue.direction == HueDirection::CLOCKWISE ? -1.0f : 1.0f;
    // Offsets, not absolute hues: differencing base-shifted hues rounds a
    // half-turn step to either side of directed_delta's SHORTEST tie.
    float relationship[4] = {};
    uint8_t relationship_count = 3;
    switch (recipe.hue.harmony) {
    case PaletteHarmony::MONOCHROMATIC:
      relationship[0] = 0.0f;
      relationship_count = 1;
      break;
    case PaletteHarmony::ANALOGOUS:
      relationship[0] = -orientation * spread;
      relationship[1] = 0.0f;
      relationship[2] = orientation * spread;
      break;
    case PaletteHarmony::ACCENTED_ANALOGOUS:
      relationship[0] = -orientation * spread;
      relationship[1] = 0.0f;
      relationship[2] = orientation * spread;
      relationship[3] = orientation * 0.5f;
      relationship_count = 4;
      break;
    case PaletteHarmony::COMPLEMENTARY:
      relationship[0] = 0.0f;
      relationship[1] = orientation * 0.5f;
      relationship_count = 2;
      break;
    case PaletteHarmony::SPLIT_COMPLEMENTARY:
      relationship[0] = 0.0f;
      relationship[1] = orientation * (0.5f - spread);
      relationship[2] = orientation * (0.5f + spread);
      break;
    case PaletteHarmony::TRIADIC:
      relationship[0] = 0.0f;
      relationship[1] = orientation / 3.0f;
      relationship[2] = orientation * 2.0f / 3.0f;
      break;
    case PaletteHarmony::TETRADIC:
      relationship[0] = 0.0f;
      relationship[1] = orientation * spread;
      relationship[2] = orientation * 0.5f;
      relationship[3] = orientation * (0.5f + spread);
      relationship_count = 4;
      break;
    case PaletteHarmony::SQUARE:
      relationship[0] = 0.0f;
      relationship[1] = orientation * 0.25f;
      relationship[2] = orientation * 0.5f;
      relationship[3] = orientation * 0.75f;
      relationship_count = 4;
      break;
    }

    for (int i = 1; i < relationship_count; ++i)
      relationship[i] = relationship[i - 1] +
                        directed_delta(relationship[i] - relationship[i - 1],
                                       recipe.hue.direction);
    for (int i = 0; i < relationship_count; ++i)
      relationship[i] += base;

    const uint8_t key_count = control_key_count(recipe);
    for (int i = 0; i < key_count; ++i) {
      const float position = i / float(key_count - 1);
      const float scaled = position * (relationship_count - 1);
      const int left = std::min(static_cast<int>(scaled),
                                static_cast<int>(relationship_count) - 1);
      const int right =
          std::min(left + 1, static_cast<int>(relationship_count) - 1);
      const float amount = scaled - left;
      hues[i] = relationship[left] +
                (relationship[right] - relationship[left]) * amount;
    }
  }

  HS_COLD_MEMBER static void resolve_hues(const PaletteRecipe &recipe,
                                          uint8_t key_count,
                                          float hues[PALETTE_MAX_KEYS],
                                          float &closing_hue) {
    if (recipe.hue.mode == HueMode::HARMONY) {
      resolve_harmony(recipe, hues);
    } else if (recipe.hue.mode == HueMode::SWEEP) {
      float sweep = recipe.hue.sweep_turns;
      if (recipe.hue.direction == HueDirection::CLOCKWISE)
        sweep = -fabsf(sweep);
      else if (recipe.hue.direction == HueDirection::COUNTERCLOCKWISE)
        sweep = fabsf(sweep);
      if (recipe.domain == PaletteDomain::LOOP) {
        for (int i = 0; i < key_count; ++i)
          hues[i] =
              recipe.hue.base_turns + sweep * i / static_cast<float>(key_count);
        closing_hue = hues[0] + sweep;
        return;
      }
      for (int i = 0; i < key_count; ++i)
        hues[i] = recipe.hue.base_turns +
                  sweep * i / static_cast<float>(key_count - 1);
    } else {
      for (int i = 0; i < key_count; ++i)
        hues[i] = recipe.hue.custom_turns[i];
    }

    const float first_raw = recipe.hue.mode == HueMode::CUSTOM
                                ? recipe.hue.custom_turns[0]
                                : hues[0];
    const float last_raw = recipe.hue.mode == HueMode::CUSTOM
                               ? recipe.hue.custom_turns[key_count - 1]
                               : hues[key_count - 1];
    closing_hue = hues[key_count - 1] +
                  directed_delta(first_raw - last_raw, recipe.hue.direction);
  }

  HS_COLD_MEMBER void initialize(const PaletteRecipe &recipe) {
    float lightness[PALETTE_MAX_KEYS];
    float chroma[PALETTE_MAX_KEYS];
    float hues[PALETTE_MAX_KEYS];
    key_count = control_key_count(recipe);
    lightness_axis = axis_state(recipe.lightness);
    chroma_axis = axis_state(recipe.chroma);
    resolve_axis(recipe.lightness, key_count, lightness);
    resolve_axis(recipe.chroma, key_count, chroma);
    resolve_hues(recipe, key_count, hues, closing_hue);

    for (int i = 0; i < key_count; ++i)
      keys[i] = {lightness[i], chroma[i], hues[i] * 2.0f * PI_F};
    for (int i = key_count; i < PALETTE_MAX_KEYS; ++i)
      keys[i] = {};
    closing_hue *= 2.0f * PI_F;
    input_offset = recipe.input.offset;
    input_span = recipe.input.span;
    domain = recipe.domain;
    easing = recipe.easing;
    color_path = recipe.color_path;
    complementary_harmony = recipe.hue.mode == HueMode::HARMONY &&
                            recipe.hue.harmony == PaletteHarmony::COMPLEMENTARY;
    chroma_basis = recipe.chroma.basis;
    headroom = recipe.chroma.headroom;
    hue_torsion = recipe.hue_torsion;
    falloff_start = recipe.falloff_start;
  }

  void assign(const Snapshot &snapshot) {
    key_count = snapshot.key_count;
    lightness_axis = {snapshot.lightness_low, snapshot.lightness_high,
                      snapshot.lightness_curve};
    chroma_axis = {snapshot.chroma_low, snapshot.chroma_high,
                   snapshot.chroma_curve};
    for (int i = 0; i < key_count; ++i)
      keys[i] = decode_snapshot_key(snapshot.keys[i]);
    for (int i = key_count; i < PALETTE_MAX_KEYS; ++i)
      keys[i] = {};
  }

  static ControlKey sample_snapshot(const Snapshot &snapshot, float position) {
    const float scaled = position * (snapshot.key_count - 1);
    const int left_index = std::min(static_cast<int>(scaled),
                                    static_cast<int>(snapshot.key_count) - 1);
    const int right_index =
        std::min(left_index + 1, static_cast<int>(snapshot.key_count) - 1);
    const float amount = scaled - left_index;
    const ControlKey left = snapshot_key(snapshot, left_index);
    const ControlKey right = snapshot_key(snapshot, right_index);
    return {left.L + (right.L - left.L) * amount,
            left.chroma + (right.chroma - left.chroma) * amount,
            left.h + (right.h - left.h) * amount};
  }

  static Snapshot::Key encode_snapshot_key(const ControlKey &key) {
    Snapshot::Key result{};
    const uint32_t lightness = static_cast<uint32_t>(
        roundf(hs::clamp(key.L, 0.0f, 1.0f) * SNAPSHOT_AXIS_STEPS));
    const uint32_t chroma = static_cast<uint32_t>(
        roundf(hs::clamp(key.chroma, 0.0f, 1.0f) * SNAPSHOT_AXIS_STEPS));
    const uint32_t axes = lightness | (chroma << 12);
    result.bytes[0] = static_cast<uint8_t>(axes);
    result.bytes[1] = static_cast<uint8_t>(axes >> 8);
    result.bytes[2] = static_cast<uint8_t>(axes >> 16);
    uint32_t hue_bits;
    std::memcpy(&hue_bits, &key.h, sizeof(hue_bits));
    for (int i = 0; i < 4; ++i)
      result.bytes[i + 3] = static_cast<uint8_t>(hue_bits >> (8 * i));
    return result;
  }

  static ControlKey decode_snapshot_key(const Snapshot::Key &key) {
    const uint32_t axes = uint32_t{key.bytes[0]} |
                          (uint32_t{key.bytes[1]} << 8) |
                          (uint32_t{key.bytes[2]} << 16);
    const float lightness = (axes & 0xFFFu) * (1.0f / SNAPSHOT_AXIS_STEPS);
    const float chroma = ((axes >> 12) & 0xFFFu) * (1.0f / SNAPSHOT_AXIS_STEPS);
    uint32_t hue_bits = 0;
    for (int i = 0; i < 4; ++i)
      hue_bits |= uint32_t{key.bytes[i + 3]} << (8 * i);
    float hue;
    std::memcpy(&hue, &hue_bits, sizeof(hue));
    return {lightness, chroma, hue};
  }

  bool is_chromatic(const ControlKey &key) const {
    return chroma_basis == ChromaBasis::LOCAL_GAMUT
               ? key.q > 0.0f
               : key.C >= OKLCH_ACHROMATIC_C;
  }

  float apply_easing(float progress) const {
    switch (easing) {
    case SegmentEase::LINEAR:
      return progress;
    case SegmentEase::COSINE:
      return 0.5f - 0.5f * fast_cosf(PI_F * progress);
    case SegmentEase::SMOOTHSTEP:
      return progress * progress * (3.0f - 2.0f * progress);
    }
    return progress;
  }

  float window_position(float position) const {
    return input_offset + hs::clamp(position, 0.0f, 1.0f) * input_span;
  }

  HS_COLD_MEMBER Segment select_base_segment(float position) const {
    position = hs::clamp(position, 0.0f, 1.0f);
    const float scaled = position * (key_count - 1);
    const int left =
        std::min(static_cast<int>(scaled), static_cast<int>(key_count) - 2);
    const int right = left + 1;
    const float progress = apply_easing(scaled - left);
    return {keys[left], keys[right], sample_envelope(position),
            1.0f,       progress,    position};
  }

  HS_COLD_MEMBER float sample_hue(float position) const {
    const Segment segment = select_base_segment(position);
    return segment.left.h +
           (segment.right.h - segment.left.h) * segment.progress;
  }

  HS_COLD_MEMBER Segment select_segment(float t) const {
    t = hs::clamp(t, 0.0f, 1.0f);

    if (domain == PaletteDomain::LOOP) {
      const float main_end = (key_count - 1.0f) / key_count;
      const float relationship = window_position(1.0f - fabsf(2.0f * t - 1.0f));
      if (t < main_end) {
        Segment segment = select_base_segment(window_position(t / main_end));
        segment.relationship_position = relationship;
        return segment;
      }

      const float progress = apply_easing((t - main_end) / (1.0f - main_end));
      const float start_position = window_position(0.0f);
      const float end_position = window_position(1.0f);
      const ControlKey start = sample_envelope(start_position);
      const ControlKey end = sample_envelope(end_position);
      const ControlKey envelope{
          end.L + (start.L - end.L) * progress,
          end.chroma + (start.chroma - end.chroma) * progress, 0.0f};
      const float winding = closing_hue - keys[0].h;
      const ControlKey left{0.0f, 0.0f, sample_hue(end_position)};
      const ControlKey right{0.0f, 0.0f, sample_hue(start_position) + winding};
      return {left, right, envelope, 1.0f, progress, relationship};
    }

    float position = t;
    float visibility = 1.0f;
    switch (domain) {
    case PaletteDomain::STRAIGHT:
      break;
    case PaletteDomain::MIRROR:
      position = 2.0f * std::min(t, 1.0f - t);
      break;
    case PaletteDomain::VIGNETTE:
      if (t < 0.1f) {
        position = 0.0f;
        visibility = apply_easing(t * 10.0f);
      } else if (t > 0.9f) {
        position = 1.0f;
        visibility = 1.0f - apply_easing((t - 0.9f) * 10.0f);
      } else {
        position = (t - 0.1f) * 1.25f;
      }
      break;
    case PaletteDomain::FALLOFF:
      if (t <= 2.0f / 3.0f) {
        position = t * 1.5f;
      } else {
        position = 1.0f;
        if (t < falloff_start)
          visibility = 1.0f - apply_easing((t - 2.0f / 3.0f) /
                                           (falloff_start - 2.0f / 3.0f));
        else
          visibility = 0.0f;
      }
      break;
    case PaletteDomain::LOOP:
      break;
    }

    Segment segment = select_base_segment(window_position(position));
    segment.visibility = visibility;
    return segment;
  }

  HS_COLD_MEMBER ControlKey sample_envelope(float position) const {
    position = hs::clamp(position, 0.0f, 1.0f);
    const float scaled = position * (key_count - 1);
    const int left =
        std::min(static_cast<int>(scaled), static_cast<int>(key_count) - 1);
    const int right = std::min(left + 1, static_cast<int>(key_count) - 1);
    const float progress = apply_easing(scaled - floorf(scaled));
    const float lightness =
        lightness_axis.curve == AxisCurve::CUSTOM
            ? keys[left].L + (keys[right].L - keys[left].L) * progress
            : evaluate_axis(lightness_axis.low, lightness_axis.high,
                            lightness_axis.curve, position);
    const float chroma =
        chroma_axis.curve == AxisCurve::CUSTOM
            ? keys[left].chroma +
                  (keys[right].chroma - keys[left].chroma) * progress
            : evaluate_axis(chroma_axis.low, chroma_axis.high,
                            chroma_axis.curve, position);
    return {lightness, chroma, 0.0f};
  }

  HS_COLD_MEMBER float realized_chroma(const ControlKey &key,
                                       float h_final) const {
    if (chroma_basis == ChromaBasis::LOCAL_GAMUT)
      return std::min(key.q, headroom) *
             gamut_continuous_chroma(key.L, h_final);
    return std::max(0.0f, key.C);
  }

  HS_COLD_MEMBER OKLCH resolve_key(const ControlKey &key,
                                   bool chromatic) const {
    const float h_path = chromatic ? key.h : 0.0f;
    const float h_final = h_path + hue_torsion * (key.L - 0.5f);
    return {key.L, chromatic ? realized_chroma(key, h_final) : 0.0f, h_final};
  }

  HS_COLD_MEMBER PathEvaluation
  evaluate_cartesian_path(const ControlKey &left_key, bool left_chromatic,
                          const ControlKey &right_key, bool right_chromatic,
                          float progress, float control) const {
    const OKLab left = oklch_to_oklab(resolve_key(left_key, left_chromatic));
    const OKLab right = oklch_to_oklab(resolve_key(right_key, right_chromatic));
    const OKLab lab{left.L + (right.L - left.L) * progress,
                    left.a + (right.a - left.a) * progress,
                    left.b + (right.b - left.b) * progress};
    const OKLCH lch = oklab_to_oklch(lab);
    const float h_path = lch.C >= OKLCH_ACHROMATIC_C ? lch.h : 0.0f;
    const float boundary = gamut_continuous_chroma(lch.L, h_path);
    const float q = boundary > 0.0f ? lch.C / boundary : control;
    return {lab, q, h_path, h_path, boundary};
  }

  HS_COLD_MEMBER PathEvaluation
  evaluate_arc_path(float L, float control, float left_h, bool left_chromatic,
                    float right_h, bool right_chromatic, float progress) const {
    float h_path;
    if (left_chromatic && right_chromatic)
      h_path = left_h + (right_h - left_h) * progress;
    else if (left_chromatic)
      h_path = left_h;
    else if (right_chromatic)
      h_path = right_h;
    else
      h_path = 0.0f;

    const float h_final = h_path + hue_torsion * (L - 0.5f);
    const float boundary = gamut_continuous_chroma(L, h_final);
    const float C = chroma_basis == ChromaBasis::LOCAL_GAMUT
                        ? std::min(control, headroom) * boundary
                        : std::max(0.0f, control);
    const float q = chroma_basis == ChromaBasis::LOCAL_GAMUT
                        ? control
                        : (boundary > 0.0f ? C / boundary : 0.0f);
    return {oklch_to_oklab({L, C, h_final}), q, h_path, h_final, boundary};
  }

  HS_COLD_MEMBER Evaluated finish_evaluation(const PathEvaluation &path) const {
    const OKLCH lch = oklab_to_oklch(path.lab);
    float r, g, blue;
    oklab_to_linear_rgb(path.lab, r, g, blue);
    return {path.lab,
            lch.C,
            path.q,
            path.C_max,
            path.h_path,
            path.h_final,
            !linear_rgb_in_gamut(r, g, blue)};
  }

  HS_COLD_MEMBER PathEvaluation evaluate_path(float t) const {
    const Segment segment = select_segment(t);
    const ControlKey envelope = segment.envelope;
    const float progress = segment.progress;
    const float visibility = segment.visibility;
    const float L = envelope.L * visibility;
    const float control = envelope.chroma * visibility;
    const bool chromatic = chroma_basis == ChromaBasis::LOCAL_GAMUT
                               ? control > 0.0f
                               : control >= OKLCH_ACHROMATIC_C;

    PathEvaluation path;
    if (complementary_harmony) {
      const float relationship_position =
          apply_easing(segment.relationship_position);

      if (color_path == ColorPath::OKLAB_CARTESIAN) {
        const ControlKey first{L, control, keys[0].h};
        const ControlKey opposite{L, control, keys[key_count - 1].h};
        path = evaluate_cartesian_path(first, chromatic, opposite, chromatic,
                                       relationship_position, control);
      } else {
        path = evaluate_arc_path(L, control, keys[0].h, chromatic,
                                 keys[key_count - 1].h, chromatic,
                                 relationship_position);
      }
    } else if (color_path == ColorPath::OKLAB_CARTESIAN) {
      const ControlKey left_key{L, control, segment.left.h};
      const ControlKey right_key{L, control, segment.right.h};
      path = evaluate_cartesian_path(left_key, chromatic, right_key, chromatic,
                                     progress, control);
    } else {
      path = evaluate_arc_path(L, control, segment.left.h, chromatic,
                               segment.right.h, chromatic, progress);
    }
    return path;
  }

  std::array<ControlKey, PALETTE_MAX_KEYS> keys{};
  uint8_t key_count = 3;
  AxisState lightness_axis{0.62f, 0.62f, AxisCurve::CONSTANT};
  AxisState chroma_axis{0.62f, 0.62f, AxisCurve::CONSTANT};
  float closing_hue = 0.0f;
  float headroom = 0.94f;
  float hue_torsion = 0.0f;
  float falloff_start = 0.90f;
  float input_offset = 0.0f;
  float input_span = 1.0f;
  PaletteDomain domain = PaletteDomain::STRAIGHT;
  SegmentEase easing = SegmentEase::COSINE;
  ColorPath color_path = ColorPath::OKLCH_ARC;
  ChromaBasis chroma_basis = ChromaBasis::LOCAL_GAMUT;
  bool complementary_harmony = false;
};

static_assert(sizeof(GenerativePalette) <= 160);

namespace PaletteRecipes {

inline float hue_turns(uint32_t hue) {
  return static_cast<float>(hue & 0xFFu) * (1.0f / 256.0f);
}

inline PaletteRecipe harmony(PaletteDomain domain, PaletteHarmony harmony,
                             float base_turns = 0.0f) {
  PaletteRecipe recipe;
  recipe.domain = domain;
  recipe.hue.harmony = harmony;
  recipe.hue.base_turns = base_turns;
  return recipe;
}

inline PaletteRecipe balanced_analogous(float base_turns = 0.0f) {
  return harmony(PaletteDomain::STRAIGHT, PaletteHarmony::ANALOGOUS,
                 base_turns);
}

inline PaletteRecipe profile(PaletteDomain domain, PaletteHarmony harmony,
                             AxisCurve lightness_curve, float base_turns,
                             float chroma = 0.62f) {
  PaletteRecipe recipe = PaletteRecipes::harmony(domain, harmony, base_turns);
  recipe.lightness.curve = lightness_curve;
  if (lightness_curve != AxisCurve::CONSTANT) {
    recipe.lightness.center = 0.52f;
    recipe.lightness.range = 0.72f;
  }
  recipe.chroma.center = chroma;
  return recipe;
}

HS_FLASH_MEMBER inline PaletteRecipe random_profile(PaletteDomain domain,
                                                    PaletteHarmony harmony,
                                                    AxisCurve lightness_curve,
                                                    float chroma = 0.62f) {
  return profile(domain, harmony, lightness_curve,
                 hue_turns(static_cast<uint32_t>(hs::rand_int(0, 256))),
                 chroma);
}

inline PaletteRecipe from_oklch_keys(PaletteDomain domain, OKLCH a, OKLCH b,
                                     OKLCH c) {
  b.h = a.h + wrap_angle_pi(b.h - a.h);
  c.h = b.h + wrap_angle_pi(c.h - b.h);
  PaletteRecipe recipe;
  recipe.domain = domain;
  recipe.hue.mode = HueMode::CUSTOM;
  recipe.hue.custom_turns[0] = a.h / (2.0f * PI_F);
  recipe.hue.custom_turns[1] = b.h / (2.0f * PI_F);
  recipe.hue.custom_turns[2] = c.h / (2.0f * PI_F);
  recipe.lightness.curve = AxisCurve::CUSTOM;
  recipe.lightness.custom[0] = a.L;
  recipe.lightness.custom[1] = b.L;
  recipe.lightness.custom[2] = c.L;
  recipe.chroma.basis = ChromaBasis::ABSOLUTE;
  recipe.chroma.curve = AxisCurve::CUSTOM;
  recipe.chroma.custom[0] = a.C;
  recipe.chroma.custom[1] = b.C;
  recipe.chroma.custom[2] = c.C;
  return recipe;
}

inline PaletteRecipe from_colors(PaletteDomain domain, const CPixel &a,
                                 const CPixel &b, const CPixel &c) {
  return from_oklch_keys(domain, pixel_to_oklch(Pixel(a)),
                         pixel_to_oklch(Pixel(b)), pixel_to_oklch(Pixel(c)));
}

inline PaletteRecipe isolight_spectral_loop(float base_turns = 0.0f) {
  PaletteRecipe recipe;
  recipe.domain = PaletteDomain::LOOP;
  recipe.hue.mode = HueMode::SWEEP;
  recipe.hue.base_turns = base_turns;
  recipe.hue.sweep_turns = 1.0f;
  recipe.lightness.center = 0.62f;
  recipe.chroma.center = 0.72f;
  return recipe;
}

inline PaletteRecipe tonal_monochrome(float base_turns = 0.0f) {
  PaletteRecipe recipe = harmony(PaletteDomain::STRAIGHT,
                                 PaletteHarmony::MONOCHROMATIC, base_turns);
  recipe.lightness.curve = AxisCurve::ASCENDING;
  recipe.lightness.center = 0.52f;
  recipe.lightness.range = 0.72f;
  recipe.chroma.curve = AxisCurve::BELL;
  recipe.chroma.center = 0.52f;
  recipe.chroma.range = 0.42f;
  return recipe;
}

} // namespace PaletteRecipes
