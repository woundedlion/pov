/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file palette_bindings.h
 * @brief Versioned GenerativePalette recipe bridge.
 *
 * Decodes a JS recipe object into a PaletteRecipe, compiles it through
 * GenerativePalette::try_compile() and bakes the 256-entry LUT (plus, for
 * inspectV4, per-sample diagnostics) into the module-global buffers below, which
 * cross back as typed memory views. Included only by targets/wasm/wasm.cpp.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <emscripten/bind.h>
#include <limits>
#include <type_traits>

#include "core/color/color.h"
#include "core/color/effect_palette_recipes.h"
#include "core/engine/platform.h"

using namespace emscripten;

static constexpr size_t PALETTE_LUT_BYTES = 256 * 3;
static constexpr size_t PALETTE_DIAGNOSTIC_FLOATS = 256 * 6;
// Bake targets for every PaletteOps instance and call, handed to JS as views
// over WASM linear memory rather than copies. See the memory-view contract on
// compile().
static uint8_t palette_lut[PALETTE_LUT_BYTES];
static float palette_diagnostics[PALETTE_DIAGNOSTIC_FLOATS];
static uint8_t palette_fallback[256];

/**
 * @brief JS-facing V4 palette-recipe compiler.
 * @details Carries no per-instance state; the bake targets are module-global.
 */
struct PaletteOps {
  /**
   * @brief Compiles a recipe and bakes its LUT.
   * @param input JS recipe object in the schema effectPresetsV4() emits.
   * @return {status} alone when the recipe is rejected, else {status,
   *         canonicalRecipe, lut}. See compile() for the memory-view contract
   *         the lut view rides on.
   */
  val compileAndBakeV4(const val &input) { return compile(input, false); }

  /**
   * @brief Compiles a recipe and bakes per-sample diagnostics with its LUT.
   * @param input JS recipe object in the schema effectPresetsV4() emits.
   * @return compileAndBakeV4()'s object plus diagnostics (256 x {L, C, q, C_max,
   *         h_path, h_final}) and fallback (one gamut-mapping flag per sample),
   *         under the same memory-view contract.
   */
  val inspectV4(const val &input) { return compile(input, true); }

  /**
   * @brief Publishes every effect-owned palette recipe.
   * @return JS array of {name, randomHue, recipe}, each recipe in the schema
   *         compileAndBakeV4() and inspectV4() accept.
   */
  val effectPresetsV4() {
    val output = val::array();
    int index = 0;
    for (const auto &preset : EffectPaletteRecipes::presets()) {
      val entry = val::object();
      entry.set("name", std::string(preset.name));
      entry.set("randomHue", preset.random_hue);
      entry.set("recipe", encode_recipe(preset.recipe));
      output.set(index++, entry);
    }
    return output;
  }

private:
  /**
   * @brief Decodes one enum field, rejecting a non-numeric field as
   *        INVALID_SCHEMA and anything the underlying type cannot hold as
   *        INVALID_ENUM. Each enum's own bound belongs to
   *        GenerativePalette::try_compile().
   * @details A missing or misspelled key reads as undefined, which converts to
   *          enumerator 0 rather than failing, so the numeric test comes first.
   */
  template <typename Enum>
  static bool decode_enum(const val &object, const char *name, Enum &output,
                          PaletteRecipeField field,
                          PaletteCompileStatus &status) {
    using Underlying = std::underlying_type_t<Enum>;
    static_assert(std::is_unsigned<Underlying>::value,
                  "recipe enums must have an unsigned underlying type");
    constexpr int LAST =
        static_cast<int>(std::numeric_limits<Underlying>::max());
    const val field_value = object[name];
    if (!field_value.isNumber()) {
      status.code = PaletteCompileCode::INVALID_SCHEMA;
      status.field = field;
      return false;
    }
    const int value = field_value.as<int>();
    if (value < 0 || value > LAST) {
      status.code = PaletteCompileCode::INVALID_ENUM;
      status.field = field;
      return false;
    }
    output = static_cast<Enum>(value);
    return true;
  }

  /**
   * @brief Decodes one PALETTE_MAX_KEYS custom-key array.
   * @param input JS array; a missing element decodes to NaN, which
   *        GenerativePalette::try_compile() rejects as NON_FINITE.
   * @param output Key array to fill.
   */
  static void decode_key_values(const val &input,
                                std::array<float, PALETTE_MAX_KEYS> &output) {
    for (int i = 0; i < PALETTE_MAX_KEYS; ++i)
      output[i] = input[i].as<float>();
  }

  /**
   * @brief Rejects a nested recipe block the caller left out.
   * @param block Value read off the recipe object.
   * @param field Field id reported for a missing block.
   * @param status Status written on rejection.
   * @return True when @p block can be indexed.
   * @details Indexing undefined or null throws a JS TypeError out through the
   *          binding instead of returning a status; every other type yields
   *          undefined leaves, which GenerativePalette's own validation rejects.
   */
  static bool block_present(const val &block, PaletteRecipeField field,
                            PaletteCompileStatus &status) {
    if (!block.isUndefined() && !block.isNull())
      return true;
    status.code = PaletteCompileCode::INVALID_SCHEMA;
    status.field = field;
    return false;
  }

  /**
   * @brief Decodes a JS recipe object into a PaletteRecipe.
   * @param input JS recipe object.
   * @param recipe Recipe to fill; partially written when the decode fails.
   * @param status Status written on rejection.
   * @return True when every block and enum decoded.
   * @details A schemaVersion other than PaletteRecipe::SCHEMA_VERSION is passed
   *          through as 0 for try_compile() to reject. Leaf scalars are not
   *          checked here: a missing one decodes to NaN and try_compile()'s
   *          finite validation names the field.
   */
  static bool decode_recipe(const val &input, PaletteRecipe &recipe,
                            PaletteCompileStatus &status) {
    if (!block_present(input, PaletteRecipeField::NONE, status))
      return false;
    const int schema_version = input["schemaVersion"].as<int>();
    recipe.schema_version = schema_version == PaletteRecipe::SCHEMA_VERSION
                                ? PaletteRecipe::SCHEMA_VERSION
                                : 0;
    const val recipe_input = input["input"];
    if (!block_present(recipe_input, PaletteRecipeField::INPUT_OFFSET, status))
      return false;
    recipe.input.offset = recipe_input["offset"].as<float>();
    recipe.input.span = recipe_input["span"].as<float>();
    if (!decode_enum(input, "domain", recipe.domain,
                     PaletteRecipeField::PALETTE_DOMAIN, status) ||
        !decode_enum(input, "easing", recipe.easing, PaletteRecipeField::EASING,
                     status) ||
        !decode_enum(input, "colorPath", recipe.color_path,
                     PaletteRecipeField::COLOR_PATH, status))
      return false;

    const val hue = input["hue"];
    if (!block_present(hue, PaletteRecipeField::HUE_MODE, status))
      return false;
    if (!decode_enum(hue, "mode", recipe.hue.mode, PaletteRecipeField::HUE_MODE,
                     status) ||
        !decode_enum(hue, "harmony", recipe.hue.harmony,
                     PaletteRecipeField::HARMONY, status) ||
        !decode_enum(hue, "direction", recipe.hue.direction,
                     PaletteRecipeField::HUE_DIRECTION, status))
      return false;
    recipe.hue.base_turns = hue["baseTurns"].as<float>();
    recipe.hue.spread_turns = hue["spreadTurns"].as<float>();
    recipe.hue.sweep_turns = hue["sweepTurns"].as<float>();
    const val custom_turns = hue["customTurns"];
    if (!block_present(custom_turns, PaletteRecipeField::CUSTOM_TURNS_0,
                       status))
      return false;
    decode_key_values(custom_turns, recipe.hue.custom_turns);

    const val lightness = input["lightness"];
    if (!block_present(lightness, PaletteRecipeField::LIGHTNESS_CURVE, status))
      return false;
    if (!decode_enum(lightness, "curve", recipe.lightness.curve,
                     PaletteRecipeField::LIGHTNESS_CURVE, status))
      return false;
    recipe.lightness.center = lightness["center"].as<float>();
    recipe.lightness.range = lightness["range"].as<float>();
    const val lightness_custom = lightness["custom"];
    if (!block_present(lightness_custom, PaletteRecipeField::LIGHTNESS_CUSTOM_0,
                       status))
      return false;
    decode_key_values(lightness_custom, recipe.lightness.custom);

    const val chroma = input["chroma"];
    if (!block_present(chroma, PaletteRecipeField::CHROMA_CURVE, status))
      return false;
    if (!decode_enum(chroma, "curve", recipe.chroma.curve,
                     PaletteRecipeField::CHROMA_CURVE, status) ||
        !decode_enum(chroma, "basis", recipe.chroma.basis,
                     PaletteRecipeField::CHROMA_BASIS, status))
      return false;
    recipe.chroma.center = chroma["center"].as<float>();
    recipe.chroma.range = chroma["range"].as<float>();
    recipe.chroma.headroom = chroma["headroom"].as<float>();
    const val chroma_custom = chroma["custom"];
    if (!block_present(chroma_custom, PaletteRecipeField::CHROMA_CUSTOM_0,
                       status))
      return false;
    decode_key_values(chroma_custom, recipe.chroma.custom);

    recipe.hue_torsion = input["hueTorsion"].as<float>();
    recipe.falloff_start = input["falloffStart"].as<float>();
    return true;
  }

  /**
   * @brief Encodes one custom-key array as a JS array.
   * @param input Key array to encode.
   * @return JS array of PALETTE_MAX_KEYS numbers.
   */
  static val
  encode_key_values(const std::array<float, PALETTE_MAX_KEYS> &input) {
    val output = val::array();
    for (int i = 0; i < PALETTE_MAX_KEYS; ++i)
      output.set(i, input[i]);
    return output;
  }

  /**
   * @brief Encodes a PaletteRecipe as the JS object decode_recipe() accepts.
   * @param recipe Recipe to encode.
   * @return JS recipe object; enums cross as their integer values.
   */
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

  /**
   * @brief Encodes a compile status as a JS object.
   * @param status Status to encode.
   * @return {code, field, wrappedFields, clampedFields, canonicalizedFields};
   *         the three adjustment masks cross as doubles, being 64-bit bitsets
   *         keyed by PaletteRecipeField.
   */
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

  /**
   * @brief Samples the palette once into the bake buffers.
   * @param index LUT slot to write.
   * @param palette Compiled palette to sample.
   * @param t Domain position to sample at.
   * @param inspect Also record the diagnostic and gamut-fallback flag.
   */
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

  /**
   * @brief Duplicates one baked sample into another slot.
   * @param destination Slot to overwrite.
   * @param source Slot to copy from.
   * @param inspect Also copy the diagnostic and gamut-fallback flag.
   */
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

  /**
   * @brief Fills all 256 LUT slots from a compiled palette.
   * @param palette Compiled palette to sample.
   * @param inspect Also fill the diagnostic and gamut-fallback buffers.
   * @details A mirroring domain samples the first half and reflects it, a
   *          looping domain samples 255 and repeats slot 0 at 255, so both seams
   *          are exact rather than left to sampling round-off.
   */
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

  /**
   * @brief Decodes, compiles and bakes one recipe.
   * @param input JS recipe object.
   * @param inspect Emit the diagnostics and fallback views alongside the lut.
   * @return {status} alone on rejection, else {status, canonicalRecipe, lut}
   *         plus {diagnostics, fallback} when @p inspect.
   * @details WASM memory-view contract: lut, diagnostics and fallback alias
   *          module-global WASM linear memory, they are NOT copies, and two
   *          events invalidate them. The next compileAndBakeV4()/inspectV4()
   *          call on any PaletteOps instance rebakes the same buffers in place,
   *          so an outstanding view silently reports the newer palette. With
   *          ALLOW_MEMORY_GROWTH=1, any subsequent heap growth detaches the
   *          underlying ArrayBuffer and leaves the view zero-length
   *          (buffer.byteLength === 0). A caller must therefore read or copy a
   *          view before its next call into the module.
   */
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

/**
 * @brief Registers PaletteOps with embind.
 * @details Bound as instance methods: the JS palette tool and the two-repo
 *          parity tests call through a constructed PaletteOps.
 */
static void bind_palette_ops() {
  class_<PaletteOps>("PaletteOps")
      .constructor<>()
      .function("compileAndBakeV4", &PaletteOps::compileAndBakeV4)
      .function("inspectV4", &PaletteOps::inspectV4)
      .function("effectPresetsV4", &PaletteOps::effectPresetsV4);
}
