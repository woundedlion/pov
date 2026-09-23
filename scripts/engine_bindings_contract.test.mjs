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
