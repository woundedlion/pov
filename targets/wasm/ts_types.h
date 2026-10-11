/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file ts_types.h
 * @brief Typed JavaScript payloads registered for emitted Embind declarations.
 */
#pragma once
#include <emscripten/bind.h>

#define HS_TS_PAYLOAD(Name)                                                    \
  /** @brief Distinct JavaScript payload identity for TypeScript emission. */  \
  struct Name : emscripten::val {                                              \
    /** @brief Borrows the existing JavaScript payload. */                     \
    Name(const emscripten::val &value) : emscripten::val(value) {}             \
  };

HS_TS_PAYLOAD(PixelsValue)
HS_TS_PAYLOAD(ParameterValuesValue)
HS_TS_PAYLOAD(PresetIdsValue)
HS_TS_PAYLOAD(DefinitionsValue)
HS_TS_PAYLOAD(ArenaMetricsValue)
HS_TS_PAYLOAD(SizesValue)
HS_TS_PAYLOAD(ResolutionsValue)
HS_TS_PAYLOAD(MeshVerticesValue)
HS_TS_PAYLOAD(MeshFacesValue)
HS_TS_PAYLOAD(MeshClassificationValue)
HS_TS_PAYLOAD(SolidRegistryValue)
HS_TS_PAYLOAD(SolidRecipeValue)
HS_TS_PAYLOAD(MeshBoundsValue)
HS_TS_PAYLOAD(MeshArenaMetricsValue)
HS_TS_PAYLOAD(PaletteResultValue)
HS_TS_PAYLOAD(PalettePresetsValue)
HS_TS_PAYLOAD(ObjectValue)
HS_TS_PAYLOAD(SnapshotValue)
HS_TS_PAYLOAD(SnapshotInputValue)
HS_TS_PAYLOAD(ChainEntriesValue)
HS_TS_PAYLOAD(ProgramValue)
HS_TS_PAYLOAD(ParameterWritesValue)
HS_TS_PAYLOAD(ChainResultValue)
HS_TS_PAYLOAD(XYZValue)
HS_TS_PAYLOAD(RGBValue)
HS_TS_PAYLOAD(LabValue)
HS_TS_PAYLOAD(ProceduralPalettesValue)
HS_TS_PAYLOAD(EffectIdentitiesValue)

#undef HS_TS_PAYLOAD

/** @brief Registers structured payload types for TypeScript emission. */
inline void bind_ts_types() {
  emscripten::register_type<PixelsValue>("PixelsValue", R"(Uint16Array)");
  emscripten::register_type<ParameterValuesValue>("ParameterValuesValue",
                                                  R"(Float32Array)");
  emscripten::register_type<PresetIdsValue>("PresetIdsValue", R"(string[])");
  emscripten::register_type<DefinitionsValue>(
      "DefinitionsValue",
      R"(import("./engine_types.js").ParameterDefinition[])");
  emscripten::register_type<ArenaMetricsValue>(
      "ArenaMetricsValue", R"(import("./engine_types.js").ArenaMetrics)");
  emscripten::register_type<SizesValue>("SizesValue",
                                        R"(Record<string, number>)");
  emscripten::register_type<ResolutionsValue>("ResolutionsValue",
                                              R"([number, number][])");
  emscripten::register_type<MeshVerticesValue>("MeshVerticesValue",
                                               R"(Float32Array | null)");
  emscripten::register_type<MeshFacesValue>(
      "MeshFacesValue", R"(import("./engine_types.js").MeshFaces | null)");
  emscripten::register_type<MeshClassificationValue>("MeshClassificationValue",
                                                     R"(Int32Array | null)");
  emscripten::register_type<SolidRegistryValue>(
      "SolidRegistryValue",
      R"(import("./engine_types.js").SolidRegistryEntry[])");
  emscripten::register_type<SolidRecipeValue>(
      "SolidRecipeValue", R"(import("./engine_types.js").SolidRecipe | null)");
  emscripten::register_type<MeshBoundsValue>(
      "MeshBoundsValue",
      R"({max_v: number; v_name: string; max_f: number; f_name: string; max_i: number; i_name: string} | null)");
  emscripten::register_type<MeshArenaMetricsValue>(
      "MeshArenaMetricsValue",
      R"(import("./engine_types.js").MeshArenaMetrics)");
  emscripten::register_type<PaletteResultValue>(
      "PaletteResultValue",
      R"(import("./engine_types.js").PaletteCompileResult)");
  emscripten::register_type<PalettePresetsValue>(
      "PalettePresetsValue",
      R"(import("./engine_types.js").PaletteEffectPreset[])");
  emscripten::register_type<ObjectValue>("ObjectValue", R"(object)");
  emscripten::register_type<SnapshotValue>(
      "SnapshotValue", R"(import("./engine_types.js").ChainSnapshot | null)");
  emscripten::register_type<SnapshotInputValue>(
      "SnapshotInputValue", R"(import("./engine_types.js").ChainSnapshot)");
  emscripten::register_type<ChainEntriesValue>(
      "ChainEntriesValue", R"({instance: string; operator: string}[])");
  emscripten::register_type<ProgramValue>(
      "ProgramValue", R"({instance: string; operator: string}[] | null)");
  emscripten::register_type<ParameterWritesValue>(
      "ParameterWritesValue", R"({name: string; value: number}[])");
  emscripten::register_type<ChainResultValue>(
      "ChainResultValue",
      R"({status: import("./engine_types.js").EnumValue; code: string; entryIndex: number})");
  emscripten::register_type<XYZValue>("XYZValue",
                                      R"({x: number; y: number; z: number})");
  emscripten::register_type<RGBValue>("RGBValue",
                                      R"({r: number; g: number; b: number})");
  emscripten::register_type<LabValue>("LabValue",
                                      R"({L: number; a: number; b: number})");
  emscripten::register_type<ProceduralPalettesValue>(
      "ProceduralPalettesValue",
      R"({name: string; a: [number,number,number]; b: [number,number,number]; c: [number,number,number]; d: [number,number,number]}[])");
  emscripten::register_type<EffectIdentitiesValue>(
      "EffectIdentitiesValue",
      R"({id: string; display_name: string; legacy_ids: string[]}[])");
}
