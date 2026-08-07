/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file palette_bindings.h
 * @brief Versioned GenerativePalette recipe bridge.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <emscripten/bind.h>

#include "core/color/color.h"
#include "core/engine/platform.h"

using namespace emscripten;

static constexpr size_t PALETTE_LUT_BYTES = 256 * 3;
static constexpr size_t PALETTE_DIAGNOSTIC_FLOATS = 256 * 6;
static uint8_t palette_lut[PALETTE_LUT_BYTES];
static float palette_diagnostics[PALETTE_DIAGNOSTIC_FLOATS];
static uint8_t palette_fallback[256];

struct PaletteOps {
  val compileAndBakeV3(const val &input) { return compile(input, false); }

  val inspectV3(const val &input) { return compile(input, true); }

private:
  template <typename Enum>
  static bool decode_enum(const val &object, const char *name, int last,
                          Enum &output, PaletteRecipeField field,
                          PaletteCompileStatus &status) {
    const int value = object[name].as<int>();
    if (value < 0 || value > last) {
      status.code = PaletteCompileCode::INVALID_ENUM;
      status.field = field;
      return false;
    }
    output = static_cast<Enum>(value);
    return true;
  }

  static void decode_key_values(const val &input,
                                std::array<float, PALETTE_MAX_KEYS> &output) {
    for (int i = 0; i < PALETTE_MAX_KEYS; ++i)
      output[i] = input[i].as<float>();
  }

  static bool decode_recipe(const val &input, PaletteRecipe &recipe,
                            PaletteCompileStatus &status) {
    const int schema_version = input["schemaVersion"].as<int>();
    recipe.schema_version = schema_version == PaletteRecipe::SCHEMA_VERSION
                                ? PaletteRecipe::SCHEMA_VERSION
                                : 0;
    recipe.key_count = static_cast<uint8_t>(input["keyCount"].as<int>());
    const val recipe_input = input["input"];
    recipe.input.offset = recipe_input["offset"].as<float>();
    recipe.input.span = recipe_input["span"].as<float>();
    if (!decode_enum(input, "domain", static_cast<int>(PaletteDomain::LOOP),
                     recipe.domain, PaletteRecipeField::PALETTE_DOMAIN,
                     status) ||
        !decode_enum(input, "easing", static_cast<int>(SegmentEase::SMOOTHSTEP),
                     recipe.easing, PaletteRecipeField::EASING, status) ||
        !decode_enum(input, "colorPath",
                     static_cast<int>(ColorPath::OKLAB_CARTESIAN),
                     recipe.color_path, PaletteRecipeField::COLOR_PATH, status))
      return false;

    const val hue = input["hue"];
    if (!decode_enum(hue, "mode", static_cast<int>(HueMode::CUSTOM),
                     recipe.hue.mode, PaletteRecipeField::HUE_MODE, status) ||
        !decode_enum(hue, "harmony", static_cast<int>(PaletteHarmony::TRIADIC),
                     recipe.hue.harmony, PaletteRecipeField::HARMONY, status) ||
        !decode_enum(
            hue, "direction", static_cast<int>(HueDirection::COUNTERCLOCKWISE),
            recipe.hue.direction, PaletteRecipeField::HUE_DIRECTION, status))
      return false;
    recipe.hue.base_turns = hue["baseTurns"].as<float>();
    recipe.hue.spread_turns = hue["spreadTurns"].as<float>();
    recipe.hue.sweep_turns = hue["sweepTurns"].as<float>();
    decode_key_values(hue["customTurns"], recipe.hue.custom_turns);

    const val lightness = input["lightness"];
    if (!decode_enum(lightness, "curve", static_cast<int>(AxisCurve::CUSTOM),
                     recipe.lightness.curve,
                     PaletteRecipeField::LIGHTNESS_CURVE, status))
      return false;
    recipe.lightness.center = lightness["center"].as<float>();
    recipe.lightness.range = lightness["range"].as<float>();
    decode_key_values(lightness["custom"], recipe.lightness.custom);

    const val chroma = input["chroma"];
    if (!decode_enum(chroma, "curve", static_cast<int>(AxisCurve::CUSTOM),
                     recipe.chroma.curve, PaletteRecipeField::CHROMA_CURVE,
                     status) ||
        !decode_enum(chroma, "basis", static_cast<int>(ChromaBasis::ABSOLUTE),
                     recipe.chroma.basis, PaletteRecipeField::CHROMA_BASIS,
                     status))
      return false;
    recipe.chroma.center = chroma["center"].as<float>();
    recipe.chroma.range = chroma["range"].as<float>();
    recipe.chroma.headroom = chroma["headroom"].as<float>();
    decode_key_values(chroma["custom"], recipe.chroma.custom);

    recipe.hue_torsion = input["hueTorsion"].as<float>();
    recipe.falloff_start = input["falloffStart"].as<float>();
    return true;
  }

  static val
  encode_key_values(const std::array<float, PALETTE_MAX_KEYS> &input) {
    val output = val::array();
    for (int i = 0; i < PALETTE_MAX_KEYS; ++i)
      output.set(i, input[i]);
    return output;
  }

  static val encode_recipe(const PaletteRecipe &recipe) {
    val recipe_input = val::object();
    recipe_input.set("offset", recipe.input.offset);
    recipe_input.set("span", recipe.input.span);

    val hue = val::object();
    hue.set("mode", static_cast<int>(recipe.hue.mode));
    hue.set("harmony", static_cast<int>(recipe.hue.harmony));
    hue.set("direction", static_cast<int>(recipe.hue.direction));
    hue.set("baseTurns", recipe.hue.base_turns);
    hue.set("spreadTurns", recipe.hue.spread_turns);
    hue.set("sweepTurns", recipe.hue.sweep_turns);
    hue.set("customTurns", encode_key_values(recipe.hue.custom_turns));

    val lightness = val::object();
    lightness.set("curve", static_cast<int>(recipe.lightness.curve));
    lightness.set("center", recipe.lightness.center);
    lightness.set("range", recipe.lightness.range);
    lightness.set("custom", encode_key_values(recipe.lightness.custom));

    val chroma = val::object();
    chroma.set("curve", static_cast<int>(recipe.chroma.curve));
    chroma.set("basis", static_cast<int>(recipe.chroma.basis));
    chroma.set("center", recipe.chroma.center);
    chroma.set("range", recipe.chroma.range);
    chroma.set("headroom", recipe.chroma.headroom);
    chroma.set("custom", encode_key_values(recipe.chroma.custom));

    val output = val::object();
    output.set("schemaVersion", recipe.schema_version);
    output.set("keyCount", recipe.key_count);
    output.set("input", recipe_input);
    output.set("domain", static_cast<int>(recipe.domain));
    output.set("easing", static_cast<int>(recipe.easing));
    output.set("colorPath", static_cast<int>(recipe.color_path));
    output.set("hue", hue);
    output.set("lightness", lightness);
    output.set("chroma", chroma);
    output.set("hueTorsion", recipe.hue_torsion);
    output.set("falloffStart", recipe.falloff_start);
    return output;
  }

  static val encode_status(const PaletteCompileStatus &status) {
    val output = val::object();
    output.set("code", static_cast<int>(status.code));
    output.set("field", static_cast<int>(status.field));
    output.set("wrappedFields",
               static_cast<double>(status.adjustments.wrapped_fields));
    output.set("clampedFields",
               static_cast<double>(status.adjustments.clamped_fields));
    output.set("canonicalizedFields",
               static_cast<double>(status.adjustments.canonicalized_fields));
    return output;
  }

  static void store_sample(int index, const GenerativePalette &palette, float t,
                           bool inspect) {
    const CRGB color = static_cast<CRGB>(palette.get(t));
    palette_lut[index * 3] = color.r;
    palette_lut[index * 3 + 1] = color.g;
    palette_lut[index * 3 + 2] = color.b;
    if (!inspect)
      return;

    const GenerativePalette::Diagnostic diagnostic = palette.diagnose(t);
    const int offset = index * 6;
    palette_diagnostics[offset] = diagnostic.L;
    palette_diagnostics[offset + 1] = diagnostic.C;
    palette_diagnostics[offset + 2] = diagnostic.q;
    palette_diagnostics[offset + 3] = diagnostic.C_max;
    palette_diagnostics[offset + 4] = diagnostic.h_path;
    palette_diagnostics[offset + 5] = diagnostic.h_final;
    palette_fallback[index] = diagnostic.fallback_mapped ? 1 : 0;
  }

  static void copy_sample(int destination, int source, bool inspect) {
    for (int channel = 0; channel < 3; ++channel)
      palette_lut[destination * 3 + channel] =
          palette_lut[source * 3 + channel];
    if (!inspect)
      return;
    for (int field = 0; field < 6; ++field)
      palette_diagnostics[destination * 6 + field] =
          palette_diagnostics[source * 6 + field];
    palette_fallback[destination] = palette_fallback[source];
  }

  static void bake(const GenerativePalette &palette, bool inspect) {
    int sample_count = 256;
    if (palette.mirrors_domain())
      sample_count = 128;
    else if (palette.loops_domain())
      sample_count = 255;

    for (int i = 0; i < sample_count; ++i)
      store_sample(i, palette, i / 255.0f, inspect);
    if (palette.mirrors_domain()) {
      for (int i = 0; i < sample_count; ++i)
        copy_sample(255 - i, i, inspect);
    } else if (palette.loops_domain()) {
      copy_sample(255, 0, inspect);
    }
  }

  static val compile(const val &input, bool inspect) {
    PaletteRecipe recipe;
    PaletteCompileStatus status;
    val output = val::object();
    if (!decode_recipe(input, recipe, status)) {
      output.set("status", encode_status(status));
      return output;
    }

    GenerativePalette palette;
    PaletteRecipe canonical;
    if (!GenerativePalette::try_compile(recipe, palette, canonical, status)) {
      output.set("status", encode_status(status));
      return output;
    }

    bake(palette, inspect);
    output.set("status", encode_status(status));
    output.set("canonicalRecipe", encode_recipe(canonical));
    output.set("lut", val(typed_memory_view(PALETTE_LUT_BYTES, palette_lut)));
    if (inspect) {
      output.set("diagnostics", val(typed_memory_view(PALETTE_DIAGNOSTIC_FLOATS,
                                                      palette_diagnostics)));
      output.set("fallback",
                 val(typed_memory_view(size_t{256}, palette_fallback)));
    }
    return output;
  }
};

static void bind_palette_ops() {
  class_<PaletteOps>("PaletteOps")
      .constructor<>()
      .function("compileAndBakeV3", &PaletteOps::compileAndBakeV3)
      .function("inspectV3", &PaletteOps::inspectV3);
}
