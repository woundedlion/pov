import { readFileSync } from 'node:fs';
import assert from 'node:assert/strict';
import { test } from 'node:test';

const source = readFileSync(
  new URL('../targets/wasm/engine_bindings.h', import.meta.url), 'utf8');

const methodBody = (signature) => {
  const signatureAt = source.indexOf(signature);
  assert.notEqual(signatureAt, -1, `${signature} is missing`);
  const open = source.indexOf('{', signatureAt);
  let depth = 1;
  for (let i = open + 1; i < source.length; i++) {
    if (source[i] === '{') depth++;
    if (source[i] === '}' && --depth === 0) return source.slice(open + 1, i);
  }
  assert.fail(`${signature} has no closing brace`);
};

test('effect installation validates parameter capacity before success', () => {
  const install = methodBody('EffectSetResult setEffect');
  const initializedAt = install.indexOf('current_effect->init();');
  const checkedAt = install.indexOf('check_param_capacity();');
  const installedAt = install.indexOf('return EffectSetResult::INSTALLED;');

  assert.ok(initializedAt >= 0);
  assert.ok(checkedAt > initializedAt);
  assert.ok(installedAt > checkedAt);
  assert.doesNotMatch(methodBody('val getParameterDefinitions()'),
    /check_param_capacity/);
  assert.doesNotMatch(methodBody('val getParamValues()'),
    /check_param_capacity/);
});
