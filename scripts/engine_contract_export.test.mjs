import { test } from 'node:test';
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const EXPORTER = fileURLToPath(new URL('./export_engine_contract.mjs', import.meta.url));
const BINDINGS = `export interface ClassHandle {}
export interface HolosphereEngine extends ClassHandle {
  getBufferLength(): number;
  getParameterAdvisories(): {name: string; foldStrength: number}[];
}
export interface MeshOps extends ClassHandle {
  transform(): MeshOps;
}
interface EmbindModule {
  HolosphereEngine: {
    new(): HolosphereEngine;
  };
}
export type MainModule = EmbindModule;
`;

/**
 * Bound morph-sweep constants as embind reports them: float32 values widened to double.
 * @param {number} scale - Factor applied to every band, so a regenerated contract can differ.
 * @returns {string} Object-literal members for the fake module.
 */
const morphConstants = (scale) => Object.entries({
  MORPH_TRUNCATE_MIN: 0.002, MORPH_TRUNCATE_MAX: 0.995,
  MORPH_TRUNCATE_PINCH_LOWER: 0.495, MORPH_TRUNCATE_PINCH_UPPER: 0.5,
  MORPH_CHAMFER_MIN: 0.02, MORPH_CHAMFER_MAX: 0.63,
}).map(([name, value]) => `${name}: ${Math.fround(value * scale)}`).join(', ');

test('browser contracts preserve emitted signatures and derive names, values and cadence', async (t) => {
  const directory = mkdtempSync(join(tmpdir(), 'engine-contract-'));
  t.after(() => rmSync(directory, { recursive: true, force: true }));
  const modulePath = join(directory, 'engine.mjs');
  const declarations = join(directory, 'holosphere_bindings.d.ts');
  const generate = (extra = false) => {
    writeFileSync(modulePath, `
class Engine {}
Engine.prototype.getBufferLength = function() {};
Engine.prototype.getParameterAdvisories = function() { return []; };
Engine.prototype.getEffectIdentities = function() {
  return [{id: 'future-effect', display_name: 'FutureEffect'}];
};
${extra ? 'Engine.prototype.futureMethod = function() {};' : ''}
class Result { constructor(value) { this.value = value; } }
Result.APPLIED = new Result(${extra ? 73 : 0});
class MeshOps {}
MeshOps.getUnsweptOps = () => ${extra ? "['zip', 'expand']" : "['expand']"};
export default async function() {
  return {HolosphereEngine: Engine, FutureResult: Result, MeshOps, SHOW_FRAMES_PER_SECOND: ${extra ? 23 : 16},
    ${morphConstants(extra ? 0.5 : 1)}};
}
`);
    writeFileSync(declarations, BINDINGS.replace('  getBufferLength(): number;',
      '  getBufferLength(): number;' + (extra ? '\n  futureMethod(_0: number): number;' : '')));
    const result = spawnSync(process.execPath, [EXPORTER, modulePath], { encoding: 'utf8' });
    assert.equal(result.status, 0, result.stderr);
    return {
      contract: readFileSync(join(directory, 'engine_contract.mjs'), 'utf8'),
      bindings: readFileSync(declarations, 'utf8'),
      aliases: readFileSync(join(directory, 'holosphere_wasm.d.ts'), 'utf8'),
    };
  };
  const first = generate();
  assert.deepEqual(generate(), first, 'generation must be byte deterministic');
  assert.equal(first.bindings, BINDINGS.replace('transform(): MeshOps;', 'transform(): MeshOps | null;'));
  assert.match(first.contract, /FPS = 16;/);
  assert.match(first.contract, /"getParameterAdvisories"/);
  assert.match(first.bindings, /getParameterAdvisories\(\): \{name: string; foldStrength: number\}\[\];/);
  assert.match(first.contract, /"FutureResult": \{\s*"APPLIED": 0/);
  assert.match(first.contract, /"id": "future-effect"/);
  assert.doesNotMatch(first.contract, /legacy_ids/);
  const changed = generate(true);
  assert.match(changed.bindings, /futureMethod\(_0: number\): number;/);
  assert.match(changed.contract, /"futureMethod"/);
  assert.match(changed.contract, /"APPLIED": 73/);
  assert.match(changed.contract, /FPS = 23;/);
  assert.doesNotMatch(changed.bindings + changed.aliases, /\bany\b/);

  const { MORPH_SWEEP } = await import(`${pathToFileURL(join(directory, 'engine_contract.mjs'))}?changed`);
  assert.deepEqual(MORPH_SWEEP, {
    bands: {
      truncate: { t: { min: 0.001, max: 0.4975, excluded: { lower: 0.2475, upper: 0.25 } } },
      chamfer: { t: { min: 0.01, max: 0.315 } },
    },
    unswept: ['expand', 'zip'],
  }, 'bands are the shortest decimals of the bound floats; unswept ops are sorted');
  writeFileSync(modulePath, readFileSync(modulePath, 'utf8').replace(/MORPH_CHAMFER_MAX: [^,}]+,?/, ''));
  const incomplete = spawnSync(process.execPath, [EXPORTER, modulePath], { encoding: 'utf8' });
  assert.notEqual(incomplete.status, 0);
  assert.match(incomplete.stderr, /engine contract is missing MORPH_CHAMFER_MAX/);
});
