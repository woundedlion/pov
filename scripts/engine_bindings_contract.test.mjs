import { readFileSync } from 'node:fs';
import assert from 'node:assert/strict';
import { test } from 'node:test';

const source = readFileSync(
  new URL('../targets/wasm/engine_bindings.h', import.meta.url), 'utf8');

test('every exported enum binds each C++ enumerator under its own name', () => {
  const headers = [
    '../targets/wasm/engine_bindings.h',
    '../targets/wasm/palette_bindings.h',
    '../targets/wasm/mesh_ops_bindings.h',
    '../core/control/params.h',
    '../core/color/palette_recipe.h',
    '../core/render/pullback/interpreter.h',
  ].map(path => readFileSync(new URL(path, import.meta.url), 'utf8')
    .replace(/\/\*[\s\S]*?\*\/|\/\/[^\n]*/gu, ''));
  const declarations = new Map();
  for (const header of headers) {
    for (const [, name, body] of header.matchAll(
      /enum\s+class\s+(\w+)(?:\s*:\s*\w+)?\s*\{([^}]+)\}/gu,
    )) {
      declarations.set(name, body.split(',').map(entry => entry.trim())
        .filter(Boolean).map(entry => entry.split(/\s*=/u)[0]));
    }
  }
  let checked = 0;
  for (const header of headers.slice(0, 3)) {
    for (const [, type, name, body] of header.matchAll(
      /emscripten::enum_<([\w:]+)>\("(\w+)"\)([^;]*);/gu,
    )) {
      const expected = declarations.get(type.split('::').at(-1));
      assert.ok(expected, `missing declaration for ${type}`);
      const bindings = [...body.matchAll(/\.value\(\s*"(\w+)"\s*,\s*([\w:]+)\s*\)/gu)];
      assert.deepEqual(bindings.map(([, exported]) => exported).sort(),
        [...expected].sort(), name);
      for (const [, exported, value] of bindings)
        assert.equal(value, `${type}::${exported}`, name);
      ++checked;
    }
  }
  assert.equal(checked, 9);
});

test('the embind engine API preserves instance and static binding names', () => {
  const instance = [
    'setResolution', 'setEffect', 'drawFrame', 'getPixels', 'getBufferLength',
    'setDisplayCaps', 'getDisplayNorthPhi', 'getDisplaySouthPhi',
    'setParameter', 'setAnimationsPaused', 'getAnimationsPaused',
    'getPresetCount', 'getPresetIndex', 'getPresetIds', 'selectPreset',
    'selectPresetById', 'synchronizePreset', 'nextPreset', 'previousPreset',
    'setPoleLod', 'getPoleLod', 'getParameterDefinitions', 'getParamValues',
    'getParamGeneration', 'getArenaMetrics', 'getEffectSizes',
    'getEffectPresetCounts', 'getFullConfigSnapshot', 'restoreFullConfigSnapshot',
    'getFullConfigFieldDefinitions', 'setShaderChain', 'setShaderChainParameters', 'setClip', 'strobeColumns',
  ];
  const statics = ['getShaderChainCatalog', 'getSupportedResolutions', 'isLive'];
  const bindings = [...source.matchAll(
    /\.(function|class_function)\(\s*"([^"]+)"\s*,\s*&HolosphereEngine::(\w+)\)/gu,
  )];
  const allNames = [...source.matchAll(/\.(?:function|class_function)\(\s*"([^"]+)"/gu)]
    .map(([, name]) => name);
  assert.deepEqual(new Set(allNames), new Set(bindings.map(([, , name]) => name)));
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
      'getFullConfigFieldDefinitions',
    ]],
    ['HS_ENABLE_CHAIN_INTERPRETER', ['setShaderChain', 'setShaderChainParameters', 'getShaderChainCatalog']],
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
