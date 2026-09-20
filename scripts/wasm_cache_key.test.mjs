import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { mkdtempSync, readFileSync, rmSync, statSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { test } from 'node:test';
import { fileURLToPath } from 'node:url';

const SCRIPT = fileURLToPath(new URL('../cmake/wasm_cache_key.cmake', import.meta.url));
const BINARY = Buffer.from([0, 97, 115, 109, 1, 0, 0, 0]);
const digest = (binary) => createHash('sha256').update(binary).digest('hex');

function fixture(t, source) {
  const directory = mkdtempSync(join(tmpdir(), 'holosphere-wasm-cache-'));
  t.after(() => rmSync(directory, { recursive: true, force: true }));
  const glue = join(directory, 'holosphere_wasm.js');
  const wasm = join(directory, 'holosphere_wasm.wasm');
  writeFileSync(glue, source);
  writeFileSync(wasm, BINARY);
  return {
    glue,
    wasm,
    read: () => readFileSync(glue, 'utf8').replaceAll('\r\n', '\n'),
    run: () => spawnSync('cmake', [`-DGLUE=${glue}`, `-DWASM=${wasm}`, '-P', SCRIPT], {
      encoding: 'utf8',
    }),
  };
}

for (const [name, source] of [
  ['release', 'var wasmBinaryFile=new URL("holosphere_wasm.wasm",import.meta.url).href;'],
  ['debug', 'var wasmBinaryFile = new URL( "holosphere_wasm.wasm", import.meta.url).href;'],
  ['multiline', 'var wasmBinaryFile = new URL(\n  "holosphere_wasm.wasm"\n  , import.meta.url).href;'],
  ['single quotes', "var wasmBinaryFile = new URL('holosphere_wasm.wasm', import.meta.url).href;"],
]) {
  test(`cache key versions ${name} glue and preserves surrounding JavaScript`, (t) => {
    const files = fixture(t, source);
    const result = files.run();
    assert.equal(result.status, 0, result.error?.message ?? result.stdout + result.stderr);
    assert.equal(files.read(), source.replace('holosphere_wasm.wasm',
      `holosphere_wasm.wasm?v=${digest(BINARY)}`));
  });
}

test('cache key is idempotent without rewriting versioned glue', (t) => {
  const files = fixture(t, 'new URL("holosphere_wasm.wasm", import.meta.url)');
  assert.equal(files.run().status, 0);
  const contents = readFileSync(files.glue);
  const modified = statSync(files.glue, { bigint: true }).mtimeNs;
  assert.equal(files.run().status, 0);
  assert.deepEqual(readFileSync(files.glue), contents);
  assert.equal(statSync(files.glue, { bigint: true }).mtimeNs, modified);
});

test('cache key refreshes when the binary changes', (t) => {
  const files = fixture(t, 'new URL("holosphere_wasm.wasm", import.meta.url)');
  assert.equal(files.run().status, 0);
  const changed = Buffer.concat([BINARY, Buffer.from('changed')]);
  writeFileSync(files.wasm, changed);
  const result = files.run();
  assert.equal(result.status, 0, result.stdout + result.stderr);
  assert.equal(files.read(), `new URL("holosphere_wasm.wasm?v=${digest(changed)}", import.meta.url)`);
});

test('cache key rejects unexpected glue without modifying it', (t) => {
  const source = 'var unrelated = "holosphere_wasm.wasm"; new URL("another.wasm", import.meta.url);';
  const files = fixture(t, source);
  const result = files.run();
  assert.notEqual(result.status, 0);
  assert.match(result.stderr, /WASM glue has no expected binary URL/);
  assert.equal(files.read(), source);
});
