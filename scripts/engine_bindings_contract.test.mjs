import { readFileSync } from 'node:fs';
import assert from 'node:assert/strict';
import { test } from 'node:test';

const source = readFileSync(
  new URL('../targets/wasm/engine_bindings.h', import.meta.url), 'utf8');

test('every exported enum binds each C++ enumerator under its own name', () => {
  const headers = [
    '../targets/wasm/engine_bindings.h',
    '../targets/wasm/workbench_bindings.h',
    '../targets/wasm/palette_bindings.h',
    '../targets/wasm/mesh_ops_bindings.h',
    '../core/control/params.h',
    '../core/color/palette_recipe.h',
    '../core/render/pullback/interpreter.h',
    '../workbench/shader/chain_snapshot.h',
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
  for (const header of headers.slice(0, 4)) {
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
  const declared = headers.slice(0, 4)
    .reduce((count, header) => count + header.split('emscripten::enum_<').length - 1, 0);
  assert.ok(checked > 0);
  assert.equal(checked, declared);
});

test('engine bindings export their implementation names without duplicates', () => {
  const bindings = [...source.matchAll(
    /\.(function|class_function)\(\s*"([^"]+)"\s*,\s*&HolosphereEngine::(\w+)\)/gu,
  )];
  const allNames = [...source.matchAll(/\.(?:function|class_function)\(\s*"([^"]+)"/gu)]
    .map(([, name]) => name);
  assert.deepEqual(new Set(allNames), new Set(bindings.map(([, , name]) => name)));
  assert.ok(bindings.length > 0);
  assert.equal(new Set(allNames).size, bindings.length);
  for (const [, , exported, implementation] of bindings)
    assert.equal(exported, implementation);
});
