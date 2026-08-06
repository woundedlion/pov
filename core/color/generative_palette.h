/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "color/color.h"

/**
 * @brief A compiled three-key perceptual palette.
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
    ControlKey a, b, c;
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

  static_assert(sizeof(Snapshot) == 36);

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

  Snapshot snapshot() const { return {a, b, c}; }

  PaletteDomain palette_domain() const { return domain; }
  ChromaBasis palette_chroma_basis() const { return chroma_basis; }
  ColorPath palette_color_path() const { return color_path; }
  SegmentEase segment_easing() const { return easing; }
  float palette_headroom() const { return headroom; }
  float palette_hue_torsion() const { return hue_torsion; }

  bool mirrors_domain() const { return domain == PaletteDomain::MIRROR; }
  bool loops_domain() const { return domain == PaletteDomain::LOOP; }

  Diagnostic diagnose(float t) const {
    const Evaluated value = evaluate(t);
    return {value.lab.L,          value.C,      value.q,
            value.C_max,          value.h_path, value.h_final,
            value.fallback_mapped};
  }

  OKLCH resolved_oklch_key(int index) const {
    const ControlKey keys[3] = {a, b, c};
    const ControlKey &key = keys[hs::clamp(index, 0, 2)];
    const float h_final = key.h + hue_torsion * (key.L - 0.5f);
    return {key.L, realized_chroma(key, h_final), h_final};
  }

  void lerp(const GenerativePalette &from, const GenerativePalette &to,
            float amount) {
    lerp(from.snapshot(), to.snapshot(), amount);
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

    const ControlKey from_keys[3] = {from.a, from.b, from.c};
    const ControlKey to_keys[3] = {to.a, to.b, to.c};
    ControlKey *out[3] = {&a, &b, &c};
    int anchor = -1;
    for (int i = 0; i < 3; ++i) {
      if (is_chromatic(from_keys[i]) && is_chromatic(to_keys[i])) {
        anchor = i;
        break;
      }
    }

    const float anchor_delta =
        anchor >= 0 ? wrap_angle_pi(to_keys[anchor].h - from_keys[anchor].h)
                    : 0.0f;
    for (int i = 0; i < 3; ++i) {
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
      *out[i] = {from_keys[i].L + (to_keys[i].L - from_keys[i].L) * amount,
                 from_keys[i].chroma +
                     (to_keys[i].chroma - from_keys[i].chroma) * amount,
                 hue};
    }
  }

  Color4 get(float t) const override {
    const LinRGB rgb = oklab_to_linear_rgb_gamut(evaluate(t).lab);
    return Color4(Pixel(float_to_pixel16(rgb.r), float_to_pixel16(rgb.g),
                        float_to_pixel16(rgb.b)),
                  1.0f);
  }

private:
  struct Unchecked {};

  struct Segment {
    ControlKey left;
    ControlKey right;
    bool left_chromatic;
    bool right_chromatic;
    float progress;
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

  static constexpr float MAX_SWEEP_TURNS = 16.0f;
  static constexpr float MAX_CUSTOM_DELTA = 16.0f;
  static constexpr float MAX_CUSTOM_ABS_INPUT = 4096.0f;
  static constexpr float MAX_ABSOLUTE_CHROMA = 1.0f;
  static constexpr float MAX_ABS_TORSION = 4.0f * PI_F;
  static constexpr float SWEEP_INTEGER_EPS = 1e-6f;

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

  static void clamp_field(float &value, float low, float high,
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
    if (!enum_in_range(recipe.hue.harmony, PaletteHarmony::TRIADIC))
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

  static bool validate_finite(const PaletteRecipe &recipe,
                              PaletteCompileStatus &status) {
    if (!finite(recipe.hue.base_turns, PaletteRecipeField::BASE_TURNS,
                status) ||
        !finite(recipe.hue.spread_turns, PaletteRecipeField::SPREAD_TURNS,
                status) ||
        !finite(recipe.hue.sweep_turns, PaletteRecipeField::SWEEP_TURNS,
                status))
      return false;
    for (int i = 0; i < 3; ++i) {
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
    for (int i = 0; i < 3; ++i) {
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
    for (int i = 0; i < 3; ++i) {
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
                  status);
  }

  static bool canonicalize(const PaletteRecipe &input, PaletteRecipe &canonical,
                           PaletteCompileStatus &status) {
    status = {};
    if (input.schema_version != PaletteRecipe::SCHEMA_VERSION)
      return fail(status, PaletteCompileCode::INVALID_SCHEMA,
                  PaletteRecipeField::SCHEMA_VERSION);
    if (!validate_enums(input, status) || !validate_finite(input, status))
      return false;

    PaletteRecipe recipe = input;
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

    for (int i = 0; i < 3; ++i) {
      if (fabsf(recipe.hue.custom_turns[i]) > MAX_CUSTOM_ABS_INPUT)
        return fail(
            status, PaletteCompileCode::HUE_LIMIT,
            static_cast<PaletteRecipeField>(
                static_cast<uint8_t>(PaletteRecipeField::CUSTOM_TURNS_0) + i));
    }
    const float custom_offset = floorf(recipe.hue.custom_turns[0]);
    if (custom_offset != 0.0f) {
      for (int i = 0; i < 3; ++i) {
        recipe.hue.custom_turns[i] -= custom_offset;
        status.adjustments.canonicalized_fields |=
            field_bit(static_cast<PaletteRecipeField>(
                static_cast<uint8_t>(PaletteRecipeField::CUSTOM_TURNS_0) + i));
      }
    }
    if (fabsf(recipe.hue.custom_turns[1] - recipe.hue.custom_turns[0]) >
        MAX_CUSTOM_DELTA)
      return fail(status, PaletteCompileCode::HUE_LIMIT,
                  PaletteRecipeField::CUSTOM_TURNS_1);
    if (fabsf(recipe.hue.custom_turns[2] - recipe.hue.custom_turns[1]) >
        MAX_CUSTOM_DELTA)
      return fail(status, PaletteCompileCode::HUE_LIMIT,
                  PaletteRecipeField::CUSTOM_TURNS_2);

    clamp_field(recipe.lightness.center, 0.0f, 1.0f,
                PaletteRecipeField::LIGHTNESS_CENTER, status);
    clamp_field(recipe.lightness.range, 0.0f, 1.0f,
                PaletteRecipeField::LIGHTNESS_RANGE, status);
    for (int i = 0; i < 3; ++i) {
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
    for (int i = 0; i < 3; ++i) {
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

  static float directed_delta(float delta, HueDirection direction) {
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

  static void resolve_axis(const AxisControls &controls, float out[3]) {
    const float low =
        hs::clamp(controls.center - controls.range * 0.5f, 0.0f, 1.0f);
    const float high =
        hs::clamp(controls.center + controls.range * 0.5f, 0.0f, 1.0f);
    switch (controls.curve) {
    case AxisCurve::CONSTANT:
      out[0] = out[1] = out[2] = controls.center;
      break;
    case AxisCurve::ASCENDING:
      out[0] = low;
      out[1] = controls.center;
      out[2] = high;
      break;
    case AxisCurve::DESCENDING:
      out[0] = high;
      out[1] = controls.center;
      out[2] = low;
      break;
    case AxisCurve::BELL:
      out[0] = out[2] = low;
      out[1] = high;
      break;
    case AxisCurve::CUP:
      out[0] = out[2] = high;
      out[1] = low;
      break;
    case AxisCurve::CUSTOM:
      out[0] = controls.custom[0];
      out[1] = controls.custom[1];
      out[2] = controls.custom[2];
      break;
    }
  }

  static void resolve_axis(const ChromaControls &controls, float out[3]) {
    AxisControls axis;
    axis.curve = controls.curve;
    axis.center = controls.center;
    axis.range = controls.range;
    for (int i = 0; i < 3; ++i)
      axis.custom[i] = controls.custom[i];
    resolve_axis(axis, out);
  }

  static void resolve_harmony(const PaletteRecipe &recipe, float hues[3]) {
    const float base = recipe.hue.base_turns;
    const float spread = recipe.hue.spread_turns;
    float raw[3] = {};
    switch (recipe.hue.harmony) {
    case PaletteHarmony::MONOCHROMATIC:
      raw[0] = raw[1] = raw[2] = base;
      break;
    case PaletteHarmony::ANALOGOUS:
      raw[0] = base - spread;
      raw[1] = base;
      raw[2] = base + spread;
      break;
    case PaletteHarmony::ACCENTED_ANALOGOUS:
      raw[0] = base - spread;
      raw[1] = base + spread;
      raw[2] = base + 0.5f;
      break;
    case PaletteHarmony::COMPLEMENTARY:
      raw[0] = base;
      raw[1] = base + 0.5f;
      raw[2] = base;
      break;
    case PaletteHarmony::SPLIT_COMPLEMENTARY:
      raw[0] = base;
      raw[1] = base + 0.5f - spread;
      raw[2] = base + 0.5f + spread;
      break;
    case PaletteHarmony::TRIADIC:
      raw[0] = base;
      raw[1] = base + 1.0f / 3.0f;
      raw[2] = base + 2.0f / 3.0f;
      break;
    }
    hues[0] = raw[0];
    for (int i = 1; i < 3; ++i)
      hues[i] = hues[i - 1] +
                directed_delta(raw[i] - raw[i - 1], recipe.hue.direction);
  }

  static void resolve_hues(const PaletteRecipe &recipe, float hues[3],
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
        hues[0] = recipe.hue.base_turns;
        hues[1] = hues[0] + sweep / 3.0f;
        hues[2] = hues[0] + 2.0f * sweep / 3.0f;
        closing_hue = hues[0] + sweep;
        return;
      }
      hues[0] = recipe.hue.base_turns;
      hues[1] = hues[0] + sweep * 0.5f;
      hues[2] = hues[0] + sweep;
    } else {
      hues[0] = recipe.hue.custom_turns[0];
      hues[1] = recipe.hue.custom_turns[1];
      hues[2] = recipe.hue.custom_turns[2];
    }

    const float first_raw = recipe.hue.mode == HueMode::CUSTOM
                                ? recipe.hue.custom_turns[0]
                                : hues[0];
    const float last_raw = recipe.hue.mode == HueMode::CUSTOM
                               ? recipe.hue.custom_turns[2]
                               : hues[2];
    closing_hue =
        hues[2] + directed_delta(first_raw - last_raw, recipe.hue.direction);
  }

  void initialize(const PaletteRecipe &recipe) {
    float lightness[3];
    float chroma[3];
    float hues[3];
    resolve_axis(recipe.lightness, lightness);
    resolve_axis(recipe.chroma, chroma);
    resolve_hues(recipe, hues, closing_hue);

    a = {lightness[0], chroma[0], hues[0] * 2.0f * PI_F};
    b = {lightness[1], chroma[1], hues[1] * 2.0f * PI_F};
    c = {lightness[2], chroma[2], hues[2] * 2.0f * PI_F};
    closing_hue *= 2.0f * PI_F;
    domain = recipe.domain;
    easing = recipe.easing;
    color_path = recipe.color_path;
    chroma_basis = recipe.chroma.basis;
    headroom = recipe.chroma.headroom;
    hue_torsion = recipe.hue_torsion;
    falloff_start = recipe.falloff_start;
  }

  void assign(const Snapshot &snapshot) {
    a = snapshot.a;
    b = snapshot.b;
    c = snapshot.c;
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

  Segment select_segment(float t) const {
    t = hs::clamp(t, 0.0f, 1.0f);
    if (domain == PaletteDomain::LOOP && t == 1.0f)
      t = 0.0f;
    if (domain == PaletteDomain::MIRROR)
      t = std::min(t, 1.0f - t);

    ControlKey stops[5] = {};
    bool chromatic[5] = {};
    float positions[5] = {};
    int count = 0;
    const ControlKey black{0.0f, 0.0f, 0.0f};

    switch (domain) {
    case PaletteDomain::STRAIGHT:
      stops[0] = a;
      stops[1] = b;
      stops[2] = c;
      positions[0] = 0.0f;
      positions[1] = 0.5f;
      positions[2] = 1.0f;
      count = 3;
      break;
    case PaletteDomain::MIRROR:
      stops[0] = a;
      stops[1] = b;
      stops[2] = c;
      positions[0] = 0.0f;
      positions[1] = 0.25f;
      positions[2] = 0.5f;
      count = 3;
      break;
    case PaletteDomain::VIGNETTE:
      stops[0] = black;
      stops[1] = a;
      stops[2] = b;
      stops[3] = c;
      stops[4] = black;
      positions[0] = 0.0f;
      positions[1] = 0.1f;
      positions[2] = 0.5f;
      positions[3] = 0.9f;
      positions[4] = 1.0f;
      count = 5;
      break;
    case PaletteDomain::FALLOFF:
      stops[0] = a;
      stops[1] = b;
      stops[2] = c;
      stops[3] = black;
      stops[4] = black;
      positions[0] = 0.0f;
      positions[1] = 1.0f / 3.0f;
      positions[2] = 2.0f / 3.0f;
      positions[3] = falloff_start;
      positions[4] = 1.0f;
      count = 5;
      break;
    case PaletteDomain::LOOP:
      stops[0] = a;
      stops[1] = b;
      stops[2] = c;
      stops[3] = {a.L, a.chroma, closing_hue};
      positions[0] = 0.0f;
      positions[1] = 1.0f / 3.0f;
      positions[2] = 2.0f / 3.0f;
      positions[3] = 1.0f;
      count = 4;
      break;
    }

    for (int i = 0; i < count; ++i)
      chromatic[i] = is_chromatic(stops[i]);
    if (domain == PaletteDomain::VIGNETTE) {
      chromatic[0] = false;
      chromatic[4] = false;
    } else if (domain == PaletteDomain::FALLOFF) {
      chromatic[3] = false;
      chromatic[4] = false;
    }

    int index = count - 2;
    for (int i = 0; i < count - 1; ++i) {
      if (t >= positions[i] && t < positions[i + 1]) {
        index = i;
        break;
      }
    }
    const float width = positions[index + 1] - positions[index];
    const float progress =
        width > 0.0f ? apply_easing(hs::clamp((t - positions[index]) / width,
                                              0.0f, 1.0f))
                     : 0.0f;
    return {stops[index], stops[index + 1], chromatic[index],
            chromatic[index + 1], progress};
  }

  float realized_chroma(const ControlKey &key, float h_final) const {
    if (chroma_basis == ChromaBasis::LOCAL_GAMUT)
      return std::min(key.q, headroom) * gamut_max_chroma(key.L, h_final);
    return std::max(0.0f, key.C);
  }

  OKLCH resolve_key(const ControlKey &key, bool chromatic) const {
    const float h_path = chromatic ? key.h : 0.0f;
    const float h_final = h_path + hue_torsion * (key.L - 0.5f);
    return {key.L, chromatic ? realized_chroma(key, h_final) : 0.0f, h_final};
  }

  Evaluated evaluate(float t) const {
    const Segment segment = select_segment(t);
    const float progress = segment.progress;
    const float control =
        segment.left.chroma +
        (segment.right.chroma - segment.left.chroma) * progress;

    OKLab lab;
    float h_path;
    float h_final;
    float q;
    if (color_path == ColorPath::OKLAB_CARTESIAN) {
      const OKLab left =
          oklch_to_oklab(resolve_key(segment.left, segment.left_chromatic));
      const OKLab right =
          oklch_to_oklab(resolve_key(segment.right, segment.right_chromatic));
      lab = {left.L + (right.L - left.L) * progress,
             left.a + (right.a - left.a) * progress,
             left.b + (right.b - left.b) * progress};
      const OKLCH lch = oklab_to_oklch(lab);
      h_path = lch.C >= OKLCH_ACHROMATIC_C ? lch.h : 0.0f;
      h_final = h_path;
      const float boundary = gamut_max_chroma(lch.L, h_final);
      q = boundary > 0.0f ? lch.C / boundary : control;
    } else {
      const float L =
          segment.left.L + (segment.right.L - segment.left.L) * progress;
      if (segment.left_chromatic && segment.right_chromatic)
        h_path = segment.left.h + (segment.right.h - segment.left.h) * progress;
      else if (segment.left_chromatic)
        h_path = segment.left.h;
      else if (segment.right_chromatic)
        h_path = segment.right.h;
      else
        h_path = 0.0f;
      h_final = h_path + hue_torsion * (L - 0.5f);
      const float boundary = gamut_max_chroma(L, h_final);
      const float C = chroma_basis == ChromaBasis::LOCAL_GAMUT
                          ? std::min(control, headroom) * boundary
                          : std::max(0.0f, control);
      lab = oklch_to_oklab({L, C, h_final});
      q = chroma_basis == ChromaBasis::LOCAL_GAMUT
              ? control
              : (boundary > 0.0f ? C / boundary : 0.0f);
    }

    const OKLCH lch = oklab_to_oklch(lab);
    float r, g, blue;
    oklab_to_linear_rgb(lab, r, g, blue);
    const float boundary = gamut_max_chroma(lch.L, h_final);
    return {lab,
            lch.C,
            q,
            boundary,
            h_path,
            h_final,
            !linear_rgb_in_gamut(r, g, blue)};
  }

  ControlKey a{};
  ControlKey b{};
  ControlKey c{};
  float closing_hue = 0.0f;
  float headroom = 0.94f;
  float hue_torsion = 0.0f;
  float falloff_start = 0.90f;
  PaletteDomain domain = PaletteDomain::STRAIGHT;
  SegmentEase easing = SegmentEase::COSINE;
  ColorPath color_path = ColorPath::OKLCH_ARC;
  ChromaBasis chroma_basis = ChromaBasis::LOCAL_GAMUT;
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

inline PaletteRecipe random_profile(PaletteDomain domain,
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
