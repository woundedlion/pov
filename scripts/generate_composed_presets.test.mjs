import { test } from 'node:test';
import assert from 'node:assert/strict';
import { copyFile, mkdir, mkdtemp, readFile, readdir, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { SHADER_DOCUMENT_EFFECTS } from './composed_effect_roster.mjs';
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
  const ids = sources.filter(source => source.includes('DESCRIPTOR_DIGEST'))
    .map(source => /EFFECT_ID = "([^"]+)"/.exec(source)[1]).sort();
  assert.deepEqual(SHADER_DOCUMENT_EFFECTS, ids);
  assert.ok(Object.isFrozen(SHADER_DOCUMENT_EFFECTS));
});

async function fixture(t) {
  const root = await mkdtemp(join(tmpdir(), 'holosphere-composed-'));
  t.after(() => rm(root, { recursive: true, force: true }));
  for (const directory of ['effects', 'patterns', 'scripts', 'targets'])
    await mkdir(join(root, directory));
  for (const [name, file, id] of [
    ['Zulu.h', 'KaleidoscopeSmooth.h', 'kaleidoscope_smooth'],
    ['Alpha.h', 'AlienBrain.h', 'alien_brain'],
  ]) {
    await copyFile(new URL(`../effects/${file}`, import.meta.url), join(root, 'effects', name));
    await copyFile(new URL(`../patterns/${id}.shader.json`, import.meta.url), join(root, 'patterns', `${id}.shader.json`));
  }
  await writeFile(join(root, 'effects/Handwritten.h'), 'static constexpr std::string_view EFFECT_ID = "handwritten";\n');
  await writeFile(join(root, 'effects/ignored.txt'), 'DESCRIPTOR_DIGEST\n');
  return root;
}

test('discovered composed headers generate includes and stable document IDs', async (t) => {
  const root = await fixture(t);
  assert.equal(await generate({ root }), 2);
  const includes = await readFile(join(root, 'targets/composed_effect_includes.h'), 'utf8');
  assert.deepEqual([...includes.matchAll(/#include "([^"]+)"/g)].map(match => match[1]),
    ['effects/Alpha.h', 'effects/Zulu.h']);
  const roster = await import(pathToFileURL(join(root, 'scripts/composed_effect_roster.mjs')));
  assert.deepEqual(roster.SHADER_DOCUMENT_EFFECTS, ['alien-brain', 'kaleidoscope-smooth']);
  assert.ok(Object.isFrozen(roster.SHADER_DOCUMENT_EFFECTS));
  assert.equal(await generate({ root, check: true }), 2);
});

test('check reports stale and missing generated artifacts without writing them', async (t) => {
  const root = await fixture(t);
  await generate({ root });
  for (const file of ['effects/Alpha.h', 'targets/composed_effect_includes.h', 'scripts/composed_effect_roster.mjs']) {
    const path = join(root, file);
    const original = await readFile(path, 'utf8');
    const stale = file.startsWith('effects/') ? original.replace(/DESCRIPTOR_DIGEST = "[^"]+"/, 'DESCRIPTOR_DIGEST = "stale"') : 'stale\n';
    await writeFile(path, stale);
    await assert.rejects(generate({ root, check: true }), error => error.message.includes(file));
    assert.equal(await readFile(path, 'utf8'), stale);
    await writeFile(path, original);
    if (file.startsWith('effects/')) continue;
    await rm(path);
    await assert.rejects(generate({ root, check: true }), error => error.message.includes(file));
    await assert.rejects(readFile(path), { code: 'ENOENT' });
    await generate({ root });
    assert.equal(await readFile(path, 'utf8'), original);
  }
});

test('removing a composed header makes both rosters stale and regeneration removes it', async (t) => {
  const root = await fixture(t);
  await generate({ root });
  await rm(join(root, 'effects/Alpha.h'));
  await assert.rejects(generate({ root, check: true }), error =>
    error.message.includes('targets/composed_effect_includes.h') &&
    error.message.includes('scripts/composed_effect_roster.mjs'));
  assert.equal(await generate({ root }), 1);
  assert.doesNotMatch(await readFile(join(root, 'targets/composed_effect_includes.h'), 'utf8'), /Alpha/);
  assert.doesNotMatch(await readFile(join(root, 'scripts/composed_effect_roster.mjs'), 'utf8'), /alien-brain/);
  assert.equal(await generate({ root, check: true }), 1);
});

test('duplicate IDs and empty composed discovery fail generation', async (t) => {
  const root = await fixture(t);
  await copyFile(join(root, 'effects/Alpha.h'), join(root, 'effects/Duplicate.h'));
  await assert.rejects(generate({ root }), /Duplicate composed effect ID: alien-brain/);
  for (const file of ['Alpha.h', 'Duplicate.h', 'Zulu.h']) await rm(join(root, 'effects', file));
  await assert.rejects(generate({ root }), /No composed presets found/);
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
