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

/** @brief JavaScript payload type for PixelsValue. */
struct PixelsValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  PixelsValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for ParameterValuesValue. */
struct ParameterValuesValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  ParameterValuesValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for PresetIdsValue. */
struct PresetIdsValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  PresetIdsValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for DefinitionsValue. */
struct DefinitionsValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  DefinitionsValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for ArenaMetricsValue. */
struct ArenaMetricsValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  ArenaMetricsValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for SizesValue. */
struct SizesValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  SizesValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for ResolutionsValue. */
struct ResolutionsValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  ResolutionsValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for MeshVerticesValue. */
struct MeshVerticesValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  MeshVerticesValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for MeshFacesValue. */
struct MeshFacesValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  MeshFacesValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for MeshClassificationValue. */
struct MeshClassificationValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  MeshClassificationValue(const emscripten::val &value)
      : emscripten::val(value) {}
};
/** @brief JavaScript payload type for SolidRegistryValue. */
struct SolidRegistryValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  SolidRegistryValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for SolidRecipeValue. */
struct SolidRecipeValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  SolidRecipeValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for MeshBoundsValue. */
struct MeshBoundsValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  MeshBoundsValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for MeshArenaMetricsValue. */
struct MeshArenaMetricsValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  MeshArenaMetricsValue(const emscripten::val &value)
      : emscripten::val(value) {}
};
/** @brief JavaScript payload type for PaletteResultValue. */
struct PaletteResultValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  PaletteResultValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for PalettePresetsValue. */
struct PalettePresetsValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  PalettePresetsValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for ObjectValue. */
struct ObjectValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  ObjectValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for SnapshotValue. */
struct SnapshotValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  SnapshotValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for SnapshotInputValue. */
struct SnapshotInputValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  SnapshotInputValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for ChainEntriesValue. */
struct ChainEntriesValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  ChainEntriesValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for ProgramValue. */
struct ProgramValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  ProgramValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for ParameterWritesValue. */
struct ParameterWritesValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  ParameterWritesValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for ChainResultValue. */
struct ChainResultValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  ChainResultValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for XYZValue. */
struct XYZValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  XYZValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for RGBValue. */
struct RGBValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  RGBValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for LabValue. */
struct LabValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  LabValue(const emscripten::val &value) : emscripten::val(value) {}
};
/** @brief JavaScript payload type for ProceduralPalettesValue. */
struct ProceduralPalettesValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  ProceduralPalettesValue(const emscripten::val &value)
      : emscripten::val(value) {}
};
/** @brief JavaScript payload type for EffectIdentitiesValue. */
struct EffectIdentitiesValue : emscripten::val {
  /**
   * @brief Wraps the existing JavaScript payload without copying it.
   * @param value Existing JavaScript payload.
   */
  EffectIdentitiesValue(const emscripten::val &value)
      : emscripten::val(value) {}
};

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
