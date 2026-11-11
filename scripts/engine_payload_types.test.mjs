import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import {
  LOADER_TYPES,
  PayloadChecker,
  boundFunctions,
  checkEngine,
  declarations,
  parseType,
  payloadAliases,
} from './engine_payload_types.mjs';

const TYPES = `
/* block */
export interface Base {
  /** Doc. */
  name: string;
}
export interface Entry extends Base {
  readonly value: number;
  mode?: 0 | 1;
  tags: Array<{id: string; on: boolean}>;
  pair: [number, string];
  table: Record<string, number>;
  pixels: Uint8Array | null;
  nested: {kind: 'a' | "b"; list: number[]}[];
}
export type Either = {kind: 'x'; x: number} | {kind: 'y'; y: string};
`;

const ENTRY = {
  name: 'n', value: 1, tags: [{ id: 'a', on: true }], pair: [1, 's'],
  table: { a: 1 }, pixels: new Uint8Array(2), nested: [{ kind: 'b', list: [1] }],
};

/** @returns {PayloadChecker} A checker over {@link TYPES}. */
const checker = () => new PayloadChecker(declarations(TYPES));

test('a value matching its declaration reaches the declaration and its bases', () => {
  const check = checker();
  assert.deepEqual(check.check(ENTRY, parseType('Entry'), 'v'), []);
  assert.deepEqual([...check.reached].sort(), ['Base', 'Entry']);
  assert.deepEqual(check.check({ ...ENTRY, mode: 1, pixels: null }, parseType('Entry'), 'v'), []);
});

test('undeclared, missing and mistyped fields are each reported', () => {
  const check = checker();
  const { name, ...unnamed } = ENTRY;
  assert.equal(name, 'n');
  assert.deepEqual(check.check({ ...ENTRY, extra: 1 }, parseType('Entry'), 'v'),
    ['v.extra: field is not declared']);
  assert.deepEqual(check.check(unnamed, parseType('Entry'), 'v'),
    ['v.name: required field is missing']);
  assert.deepEqual(check.check({ ...ENTRY, pair: [1, 2] }, parseType('Entry'), 'v'),
    ['v.pair[1]: expected string, got number']);
  assert.deepEqual(check.check({ ...ENTRY, pixels: new Uint16Array(2) }, parseType('Entry'), 'v'),
    ['v.pixels: matches no member of the union', 'v.pixels: expected Uint8Array, got Uint16Array']);
  assert.equal(check.check({ ...ENTRY, mode: 2 }, parseType('Entry'), 'v')[0],
    'v.mode: matches no member of the union');
  assert.deepEqual([...check.reached], [], 'a failed match reaches nothing');
});

test('a union reports the closest member and reaches only the matching one', () => {
  const check = checker();
  assert.deepEqual(check.check({ kind: 'y', y: 's' }, parseType('Either'), 'v'), []);
  assert.deepEqual(check.check({ kind: 'y', y: 1 }, parseType('Either'), 'v'),
    ['v: matches no member of the union', 'v.y: expected string, got number']);
  assert.deepEqual(check.check(null, parseType('Entry | null'), 'v'), []);
  assert.equal(check.check({ ...ENTRY, value: 's' }, parseType('Entry | null'), 'v')[1],
    'v.value: expected number, got string');
});

test('an undeclared type name is an error rather than a pass', () => {
  assert.throws(() => checker().check({}, parseType('Missing'), 'v'), /Missing is not declared/);
  assert.throws(() => parseType('(a: number) => void'), /expected|unexpected|bad member/);
});

test('emitted declarations yield payload aliases and owned bound functions', () => {
  const bindings = `type EmbindString = string;
type CountsValue = Record<string, number>;
export interface Thing extends ClassHandle {
  counts(): CountsValue;
  scale(_0: number, _1: number): number;
}
interface EmbindModule {
  Thing: {
    new(): Thing;
    lookup(_0: EmbindString): CountsValue;
  };
  free(_0: number): CountsValue;
}
`;
  assert.deepEqual([...payloadAliases(bindings)], ['CountsValue']);
  assert.deepEqual(boundFunctions(bindings), [
    { owner: 'Thing', name: 'counts', params: [], returns: 'CountsValue' },
    { owner: 'Thing', name: 'scale', params: ['number', 'number'], returns: 'number' },
    { owner: 'Thing.static', name: 'lookup', params: ['EmbindString'], returns: 'CountsValue' },
    { owner: 'module', name: 'free', params: ['number'], returns: 'CountsValue' },
  ]);
});

test('checkEngine flags payload drift, unsampled functions and unreached declarations', () => {
  const bindings = `type MetricsValue = import("./engine_types.js").Metrics;
type IdentitiesValue = {id: string}[];
type RegistryValue = string[];
export interface HolosphereEngine extends ClassHandle {
  getEffectIdentities(): IdentitiesValue;
  getMetrics(): MetricsValue;
}
interface EmbindModule {
  MeshOps: {
    getRegistry(): RegistryValue;
    getNamed(_0: EmbindString): MetricsValue;
  };
}
`;
  const types = `export interface Metrics { used: number; }
export interface Orphan { x: number; }
export interface RuntimeState { dead?: boolean; }
`;
  let metrics = { used: 1, spare: 2 };
  class HolosphereEngine {
    getEffectIdentities() { return [{ id: 'A' }]; }
    getMetrics() { return metrics; }
    setEffect() {}
    drawFrame() {}
    getShaderChainBindings() { return null; }
    delete() {}
  }
  class PaletteOps {
    effectPresetsV4() { return []; }
    delete() {}
  }
  const module = {
    HolosphereEngine, PaletteOps,
    MeshOps: { getRegistry: () => [], clearToolingMemory() {} },
  };
  assert.ok(LOADER_TYPES.has('RuntimeState'));
  assert.deepEqual(checkEngine(module, bindings, types), [
    'HolosphereEngine.getMetrics().spare: field is not declared',
    'MeshOps.static.getNamed: returns a payload but no sample call checks it',
    'Metrics: declared in engine_types.d.ts but no engine payload reaches it',
    'Orphan: declared in engine_types.d.ts but no engine payload reaches it',
  ]);
  metrics = { used: 1 };
  assert.deepEqual(checkEngine(module, bindings, types).filter((line) => !line.startsWith('MeshOps')),
    ['Orphan: declared in engine_types.d.ts but no engine payload reaches it']);
});

test('the authored engine_types.d.ts parses in full', () => {
  const source = readFileSync(new URL('../targets/wasm/engine_types.d.ts', import.meta.url), 'utf8');
  const declared = declarations(source);
  const check = new PayloadChecker(declared);
  assert.ok(declared.size > 10);
  for (const name of declared.keys())
    if (!LOADER_TYPES.has(name)) assert.ok(check.resolve(name), name);
  for (const name of LOADER_TYPES) assert.ok(declared.has(name), name);
});
