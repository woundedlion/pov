import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile, readdir } from 'node:fs/promises';
import { departureLiteral, floatLiteral, generate, generatedSections, presetAssignments, updateHeader } from './generate_composed_presets.mjs';
import { compileShaderDocument } from './shader_workbench.mjs';
import { loadOperatorCatalog } from './pattern_documents.mjs';

const catalog = await loadOperatorCatalog();
const document = JSON.parse(await readFile(new URL('../patterns/kaleidoscope_smooth.shader.json', import.meta.url), 'utf8'));

test('composed preset headers reproduce their document values and digests', async () => {
  const effects = new URL('../effects/', import.meta.url);
  const headers = (await readdir(effects)).filter(name => name.endsWith('.h'));
  const sources = await Promise.all(headers.map(name => readFile(new URL(name, effects), 'utf8')));
  const count = sources.filter(source => source.includes('DESCRIPTOR_DIGEST')).length;
  assert.equal(await generate({ check: true }), count);
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

test('departure edges must target the generated successor', () => {
  const bank = structuredClone(document.preset_bank);
  bank.edges[0].to = bank.edges[0].from;
  assert.throws(() => departureLiteral(bank, bank.edges[0].from, document.descriptor.path_policies), /successor/);
});

test('departure policy kinds must be parallel', () => {
  const bank = structuredClone(document.preset_bank);
  const policies = structuredClone(document.descriptor.path_policies);
  policies[0].kind = 'SERIAL';
  assert.throws(() => departureLiteral(bank, bank.edges[0].from, policies), /policy kind/);
});

test('missing departure edges require automatic snap fallback', () => {
  const bank = structuredClone(document.preset_bank);
  bank.edges = [];
  const id = bank.choreography.generated_order[0];
  assert.throws(() => departureLiteral(bank, id), /fall back to SNAP/);
  bank.absent_edge_fallback.automatic = 'SNAP';
  assert.equal(departureLiteral(bank, id), 'Segue::Preset::Snap{}');
});


test('composed topology cannot vary between presets', () => {
  const edited = structuredClone(document);
  assert.ok(edited.preset_bank.presets.length >= 2);
  const id = 'colorize.palette-mode';
  const parameter = edited.descriptor.parameters.find((entry) => entry.id === id);
  const first = edited.preset_bank.presets[0].values[id];
  edited.preset_bank.presets[1].values[id] = parameter.domain.values.find((value) => value !== first);
  const compiled = compileShaderDocument(edited, { catalog });
  assert.equal(compiled.status, 'VALID');
  assert.throws(() => generatedSections(compiled), /topology must be uniform: colorize.palette-mode/);
});
