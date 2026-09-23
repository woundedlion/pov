import { readFileSync } from 'node:fs';
import assert from 'node:assert/strict';
import { test } from 'node:test';

const source = readFileSync(
  new URL('../targets/wasm/engine_bindings.h', import.meta.url), 'utf8');

test('the embind engine API preserves instance and static binding names', () => {
  const instance = [
    'setResolution', 'setEffect', 'drawFrame', 'getPixels', 'getBufferLength',
    'setParameter', 'setAnimationsPaused', 'getAnimationsPaused',
    'getPresetCount', 'getPresetIndex', 'getPresetIds', 'selectPreset',
    'selectPresetById', 'synchronizePreset', 'nextPreset', 'previousPreset',
    'setPoleLod', 'getPoleLod', 'getParameterDefinitions', 'getParamValues',
    'getParamGeneration', 'getArenaMetrics', 'getEffectSizes',
    'getEffectPresetCounts', 'getFullConfigSnapshot', 'restoreFullConfigSnapshot',
    'getFullConfigFieldDefinitions', 'getConfigImportNotice',
    'clearConfigImportNotice', 'setShaderChain', 'setClip', 'strobeColumns',
  ];
  const statics = ['getShaderChainCatalog', 'getSupportedResolutions', 'isLive'];
  const bindings = [...source.matchAll(
    /\.(function|class_function)\(\s*"([^"]+)"\s*,\s*&HolosphereEngine::(\w+)\)/gu,
  )];
  assert.deepEqual(bindings.filter(([, kind]) => kind === 'function')
    .map(([, , name]) => name).sort(), instance.sort());
  assert.deepEqual(bindings.filter(([, kind]) => kind === 'class_function')
    .map(([, , name]) => name).sort(), statics.sort());
  for (const [, , exported, implementation] of bindings)
    assert.equal(exported, implementation);
});

test('optional engine APIs stay inside their feature guards', () => {
  const expected = new Map([
    ['HS_ENABLE_SHADER_WORKBENCH', [
      'getFullConfigSnapshot', 'restoreFullConfigSnapshot',
      'getFullConfigFieldDefinitions', 'getConfigImportNotice',
      'clearConfigImportNotice',
    ]],
    ['HS_ENABLE_CHAIN_INTERPRETER', ['setShaderChain', 'getShaderChainCatalog']],
  ]);
  const registration = source.slice(source.indexOf('static void bind_engine()'));
  const guarded = new Map([...expected.keys()].map(flag => [flag, []]));
  let guard = null;
  for (const line of registration.split(/\r?\n/u)) {
    const directive = line.match(/^#if (HS_ENABLE_\w+)$/u);
    if (directive) guard = directive[1];
    if (/^#endif\b/u.test(line)) guard = null;
    const binding = line.match(/\.(?:function|class_function)\("([^"]+)"/u);
    if (binding && guarded.has(guard)) guarded.get(guard).push(binding[1]);
  }
  for (const [flag, names] of expected)
    assert.deepEqual(guarded.get(flag).sort(), names.sort(), flag);
});
