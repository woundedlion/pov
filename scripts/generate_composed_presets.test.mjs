import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { floatLiteral, generate, generatedSections, presetAssignments, updateHeader } from './generate_composed_presets.mjs';
import { compileShaderDocument } from './shader_workbench.mjs';
import { loadOperatorCatalog } from './pattern_documents.mjs';

const catalog = await loadOperatorCatalog();
const document = JSON.parse(await readFile(new URL('../patterns/kaleidoscope_smooth.shader.json', import.meta.url), 'utf8'));

test('composed preset headers reproduce their document values and digests', async () => {
  assert.equal(await generate({ check: true }), 18);
});

test('float literals round-trip the exact stored float without decimal noise', () => {
  for (const value of [0, -0, 0.366, 0.0269999988, 1e-40, -0.00021158854, 3.402823466e38]) {
    const literal = floatLiteral(value);
    assert.ok(Object.is(Math.fround(Number(literal.slice(0, -1))), Math.fround(value)));
  }
  assert.equal(floatLiteral(Math.fround(0.366)), '0.366f');
  assert.throws(() => floatLiteral(Infinity), /Nonfinite/);
});

test('an authored value changes the generated assignment and bank identity', () => {
  const original = generatedSections(compileShaderDocument(document, { catalog }));
  const edited = structuredClone(document);
  edited.preset_bank.presets[0].values['sample.pattern-freq'] = 3;
  const next = generatedSections(compileShaderDocument(edited, { catalog }));
  assert.notEqual(next.identity, original.identity);
  assert.match(next.params, /value.source.pattern_freq = 3.0f;/);
  assert.throws(() => updateHeader('', compileShaderDocument(document, { catalog })), /markers/);
});

test('derived and unmapped values cannot silently enter composed presets', () => {
  assert.throws(() => presetAssignments(document, { 'unknown.value': 1 }), /Unmapped/);
  assert.throws(() => presetAssignments(document, { 'warp1.lattice-period': 1 }), /Invalid derived/);
  assert.throws(() => presetAssignments(document, { 'colorize.palette-mapping': 'unknown' }), /Invalid mapping/);
});
