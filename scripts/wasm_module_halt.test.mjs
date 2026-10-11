import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runInNewContext } from 'node:vm';

const source = readFileSync(new URL('../targets/wasm/module_halt.js', import.meta.url), 'utf8');

test('abort marks the module dead before invoking the supplied callback', () => {
  const module = {
    onAbort(reason) {
      assert.equal(this, module);
      assert.equal(this.HS_MODULE_DEAD, true);
      assert.equal(reason, 'out of memory');
    },
  };
  runInNewContext(source, { Module: module });
  assert.equal(module.HS_MODULE_DEAD, undefined);
  module.onAbort('out of memory');
  assert.equal(module.HS_MODULE_DEAD, true);
});

test('abort marks the module dead without a supplied callback', () => {
  const module = {};
  runInNewContext(source, { Module: module });
  module.onAbort('out of memory');
  assert.equal(module.HS_MODULE_DEAD, true);
});

test('a throwing abort callback leaves the module marked dead', () => {
  const error = new Error('observer failed');
  const module = { onAbort() { throw error; } };
  runInNewContext(source, { Module: module });
  assert.throws(() => module.onAbort('out of memory'), (thrown) => thrown === error);
  assert.equal(module.HS_MODULE_DEAD, true);
});
