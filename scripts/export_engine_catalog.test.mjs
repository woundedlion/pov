import assert from 'node:assert/strict';
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { test } from 'node:test';

test('catalog export validates the result before replacing the output', (t) => {
  const directory = mkdtempSync(join(tmpdir(), 'catalog-export-'));
  t.after(() => rmSync(directory, { recursive: true, force: true }));
  const modulePath = join(directory, 'engine.mjs');
  const outputPath = join(directory, 'catalog.json');
  const script = fileURLToPath(new URL('./export_engine_catalog.mjs', import.meta.url));
  const valid = JSON.stringify({ operators: [{ id: 'fixture' }], carriers: ['sphere'], budgets: {} });
  for (const catalog of [valid, 'not-json', '{}']) {
    writeFileSync(modulePath, `export default async () => ({ HolosphereEngine: {
      getShaderChainCatalog: () => ${JSON.stringify(catalog)} } });`);
    writeFileSync(outputPath, 'previous');
    // NODE_OPTIONS carries the module-roster preload into this CLI-only module.
    const result = spawnSync(process.execPath, [script, modulePath, outputPath], {
      encoding: 'utf8', env: process.env,
    });
    if (catalog === valid) {
      assert.equal(result.status, 0, result.stderr);
      assert.equal(readFileSync(outputPath, 'utf8'), `${valid}\n`);
    } else {
      assert.equal(result.status, 1);
      assert.match(result.stderr, /catalog export failed:/);
      assert.equal(readFileSync(outputPath, 'utf8'), 'previous');
    }
  }
});
