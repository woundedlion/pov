import { test } from 'node:test';
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';

const EXPORTER = fileURLToPath(new URL('../scripts/export_engine_contract.mjs', import.meta.url));
const BINDINGS = `export interface ClassHandle {}
export interface HolosphereEngine extends ClassHandle {
  getBufferLength(): number;
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

test('browser contracts preserve emitted signatures and derive names, values and cadence', (t) => {
  const directory = mkdtempSync(join(tmpdir(), 'engine-contract-'));
  t.after(() => rmSync(directory, { recursive: true, force: true }));
  const modulePath = join(directory, 'engine.mjs');
  const declarations = join(directory, 'holosphere_bindings.d.ts');
  const generate = (extra = false) => {
    writeFileSync(modulePath, `
class Engine {}
Engine.prototype.getBufferLength = function() {};
Engine.prototype.getEffectIdentities = function() {
  return [{id: 'future-effect', display_name: 'FutureEffect', legacy_ids: ['OldFuture']}];
};
${extra ? 'Engine.prototype.futureMethod = function() {};' : ''}
class Result { constructor(value) { this.value = value; } }
Result.APPLIED = new Result(${extra ? 73 : 0});
export default async function() {
  return {HolosphereEngine: Engine, FutureResult: Result, SHOW_FRAMES_PER_SECOND: ${extra ? 23 : 16}};
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
  assert.match(first.contract, /"FutureResult": \{\s*"APPLIED": 0/);
  assert.match(first.contract, /"id": "future-effect"/);
  assert.match(first.contract, /"OldFuture"/);
  const changed = generate(true);
  assert.match(changed.bindings, /futureMethod\(_0: number\): number;/);
  assert.match(changed.contract, /"futureMethod"/);
  assert.match(changed.contract, /"APPLIED": 73/);
  assert.match(changed.contract, /FPS = 23;/);
  assert.doesNotMatch(changed.bindings + changed.aliases, /\bany\b/);
});
