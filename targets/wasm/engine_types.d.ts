/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

export interface EnumValue {
  readonly value: number;
}

export interface ArenaUsage {
  usage: number;
  high_water_mark: number;
  lifetime_high_water_mark: number;
  capacity: number;
}

export interface StackUsage {
  high_water_mark: number;
  init_high_water_mark: number;
  capacity: number;
}

export interface ArenaMetrics {
  scratch_arena_a: ArenaUsage;
  scratch_arena_b: ArenaUsage;
  persistent_arena: ArenaUsage;
  stack: StackUsage;
}

export interface ParameterDefinitionBase {
  /** Parameter name, as setParameter() takes it. */
  name: string;
  /** Registered pipeline folder metadata. */
  group?: {name: string; title: string; order: number};
  /**
   * True for a parameter registered as animated. An accepted setParameter() write
   * to an animated param engages the engine's animation pause.
   */
  animated: boolean;
  /** True for engine-written telemetry; setParameter() answers it READONLY. */
  readonly: boolean;
  /** True when a preset export carries this parameter. */
  preset: boolean;
  /** Actionable reason the requested value cannot currently be rendered. */
  warning?: string;
}

export interface BooleanParameterDefinition extends ParameterDefinitionBase {
  /** Current rendered value shown by the GUI. */
  value: boolean;
  /** Writable target copied when another renderer is initialized. */
  requestedValue: boolean;
  /** Last value admitted for rendering. */
  acceptedValue: boolean;
}

export interface NumericParameterDefinition extends ParameterDefinitionBase {
  /** Current rendered value shown by the GUI. */
  value: number;
  /** Writable target copied when another renderer is initialized. */
  requestedValue: number;
  /** Last value admitted for rendering. */
  acceptedValue: number;
  min: number;
  max: number;
  step?: number;
  options?: string[];
  optionValues?: number[];
  exportOptions?: string[];
}

export type ParameterDefinition =
  | BooleanParameterDefinition
  | NumericParameterDefinition;

export interface ChainSpatialWalkState {
  noiseSeed: number;
  walkTime: number;
  angularVelocity: number;
  spinPhase: number;
  position: number[];
  direction: number[];
  wander: number[];
}

export type ChainRuntimeState =
  | {instance: string; kind: 'spatial-walk-v2'; state: ChainSpatialWalkState}
  | {instance: string; kind: 'source-clock-v1'; state: {primary: number; secondary: number; angle: number}}
  | {instance: string; kind: 'noise-clock-v1'; state: {phase: number; noiseSeed: number}}
  | {instance: string; kind: 'phase-clock-v1' | 'ripple-clock-v1'; state: {phase: number}}
  | {instance: string; kind: 'affine-clock-v1'; state: {phase: number; rotation: number}}
  | {instance: string; kind: 'color-clock-v1'; state: {oscillationPhase: number; hueNoisePhase: number; hueNoiseSeed: number}}
  | {instance: string; kind: 'spherical-rings-v2'; state: {walk: ChainSpatialWalkState; phase: number}};

export interface ChainSnapshot {
  schemaVersion: number;
  chain: Array<{instance: string; operator: string}>;
  parameters: Array<{name: string; value: number}>;
  runtime?: ChainRuntimeState[];
  paletteBank?: {chroma: number; hues: number[]; cycles: Array<{frame: number; nextSequence: number; fadeActive: boolean; displayDirty: boolean}>};
  animationsPaused: boolean;
}

export interface MeshArenaMetrics {
  scratch_arena_a: ArenaUsage;
  scratch_arena_b: ArenaUsage;
  persistent_arena: ArenaUsage;
  tooling_arena: ArenaUsage;
  tooling_scratch_a: ArenaUsage;
  tooling_scratch_b: ArenaUsage;
}

export interface SolidRegistryEntry {
  /** Name fromSolidName() and getRecipe() take. */
  name: string;
  /** `Simple` for a registry seed, `Complex` for a recipe-built solid. */
  category: string;
}

export interface SolidRecipeStep {
  op: string;
  param: number;
  twist: number;
}

export interface SolidRecipe {
  seed: string;
  ops: SolidRecipeStep[];
}

export interface MeshFaces {
  indices: Uint16Array;
  counts: Uint8Array;
}

export interface PaletteCompileStatus {
  code: EnumValue;
  /** Which recipe field an error names. */
  field: EnumValue;
  /** Bitmask over the recipe fields the compiler wrapped into range. */
  wrappedFields: number;
  /** The fields it clamped, in the same bit positions. */
  clampedFields: number;
  /** The fields it rewrote into canonical form, in the same bit positions. */
  canonicalizedFields: number;
}

export interface PaletteCompileResult {
  status: PaletteCompileStatus;
  /** The recipe as the compiler normalized it. */
  canonicalRecipe?: object;
  /** 256 sRGB triples. */
  lut?: Uint8Array;
  /** Six values per entry: L, C, q, Cmax, hPath, hFinal. */
  diagnostics?: Float32Array;
  /** Per entry, non-zero where the color was mapped back into gamut. */
  fallback?: Uint8Array;
}

export interface PaletteEffectPreset {
  name: string;
  /** True where the effect randomizes the base hue at runtime. */
  randomHue: boolean;
  recipe: object;
}

/** Runtime state set by the engine's terminal trap bridge. */
export interface RuntimeState {
  HS_MODULE_DEAD?: boolean;
  onAbort?: (reason: unknown) => void;
}
/** Emscripten loader hooks supplied by the browser. */
export interface ModuleOptions {
  print?: (text: string) => void;
  printErr?: (text: string) => void;
  instantiateWasm?: (imports: WebAssembly.Imports,
    onInstance: (instance: WebAssembly.Instance, module: WebAssembly.Module) => void) => object;
}
