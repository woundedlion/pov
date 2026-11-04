import { test } from 'node:test';
import assert from 'node:assert/strict';
import { copyFile, mkdir, mkdtemp, readFile, readdir, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { SHADER_DOCUMENT_EFFECTS } from './composed_effect_roster.mjs';
import { departureLiteral, floatLiteral, generate, generatedSections, presetAssignments, topologySelections, updateHeader } from './generate_composed_presets.mjs';
import { compileShaderDocument, declarationFromCatalogField } from './shader_workbench.mjs';
import { loadOperatorCatalog } from './pattern_documents.mjs';
import { loadEffectIdentities } from './effect_roster.mjs';

const catalog = await loadOperatorCatalog();
const identities = await loadEffectIdentities();
const document = JSON.parse(await readFile(new URL('../patterns/kaleidoscope_smooth.shader.json', import.meta.url), 'utf8'));

test('composed preset headers reproduce their document values and digests', async () => {
  const effects = new URL('../effects/', import.meta.url);
  const headers = (await readdir(effects)).filter(name => name.endsWith('.h'));
  const sources = await Promise.all(headers.map(name => readFile(new URL(name, effects), 'utf8')));
  const count = sources.filter(source => source.includes('DESCRIPTOR_DIGEST')).length;
  assert.equal(await generate({ check: true }), count);
  const ids = sources.filter(source => source.includes('DESCRIPTOR_DIGEST'))
    .map(source => identities.find(row => row.name ===
      /EFFECT_IDENTITY\s*=\s*hs::effect_identities::(\w+)/u.exec(source)[1]).id).sort();
  assert.deepEqual(SHADER_DOCUMENT_EFFECTS, ids);
  assert.ok(Object.isFrozen(SHADER_DOCUMENT_EFFECTS));
});

async function fixture(t) {
  const root = await mkdtemp(join(tmpdir(), 'holosphere-composed-'));
  t.after(() => rm(root, { recursive: true, force: true }));
  for (const directory of ['effects', 'patterns', 'scripts', 'targets'])
    await mkdir(join(root, directory));
  await mkdir(join(root, 'core/platform'), { recursive: true });
  await copyFile(new URL('../core/platform/effect_identity.h', import.meta.url),
    join(root, 'core/platform/effect_identity.h'));
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
  assert.match(next.params, /value\.template get<"source">\(\)\.pattern_freq = 3\.0f;/);
  assert.throws(() => updateHeader('', compileShaderDocument(document, { catalog })), /markers/);
});

test('derived and unmapped values cannot silently enter composed presets', () => {
  assert.throws(() => presetAssignments(document, { 'unknown.value': 1 }), /Unmapped/);
  assert.throws(() => presetAssignments(document, { 'warp1.lattice-period': 1 }), /Unmapped/);
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
  const id = 'colorize.brightness-envelope';
  const parameter = edited.descriptor.parameters.find((entry) => entry.id === id);
  const first = edited.preset_bank.presets[0].values[id];
  edited.preset_bank.presets[1].values[id] = parameter.domain.values.find((value) => value !== first);
  const compiled = compileShaderDocument(edited, { catalog });
  assert.equal(compiled.status, 'VALID');
  assert.throws(() => generatedSections(compiled), /topology must be uniform: colorize.brightness-envelope/);
});
test('a mismatched document identity cannot rewrite its composed header', async (t) => {
  const root = await fixture(t);
  const path = join(root, 'patterns/alien_brain.shader.json');
  const authored = JSON.parse(await readFile(path, 'utf8'));
  authored.effect_id = 'other-effect';
  await writeFile(path, JSON.stringify(authored));
  const header = join(root, 'effects/Alpha.h');
  const original = await readFile(header, 'utf8');
  await assert.rejects(generate({ root }), /names effect "other-effect", expected "alien-brain"/);
  assert.equal(await readFile(header, 'utf8'), original);
});

test('a composed document requires authored dwell', () => {
  const authored = structuredClone(document);
  delete authored.preset_bank.choreography.dwell;
  const compiled = compileShaderDocument(authored, { catalog });
  assert.equal(compiled.status, 'VALID');
  assert.throws(() => generatedSections(compiled), /must declare choreography.dwell/);
});


test('invalid authored preset reports its document diagnostics before generation', async (t) => {
  const root = await fixture(t);
  const path = join(root, 'patterns/kaleidoscope_smooth.shader.json');
  const authored = JSON.parse(await readFile(path, 'utf8'));
  authored.preset_bank.presets[0].values['sample.pattern-freq'] = 1e20;
  await writeFile(path, JSON.stringify(authored));
  await assert.rejects(generate({ root, check: true }), error =>
    error.message.includes('patterns/kaleidoscope_smooth.shader.json') &&
    error.message.includes('VALUE_OUT_OF_RANGE'));
});

test('GridWorms keeps its authored stage order and keyed twist controls', async () => {
  const source = JSON.parse(await readFile(new URL('../patterns/grid_worms.shader.json', import.meta.url), 'utf8'));
  const compiled = compileShaderDocument(source, { catalog });
  assert.equal(compiled.status, 'VALID');
  assert.equal(source.document_id, 'grid-worms-v1');
  assert.equal(source.effect_id, 'grid-worms');
  assert.deepEqual(source.descriptor.chain.map(({ operator }) => operator), [
    'sphere.lens.twist.v2', 'sphere.rotate.v2', 'project.gnomonic.v2',
    'sample.lattice.v2', 'colorize.generated-palette.v4', 'color.modulate.v1',
  ]);
  const values = source.preset_bank.presets[0].values;
  const assignments = presetAssignments(source, values);
  assert.equal(assignments.get('lens1.twist_angle'), '-10.259509f');
  assert.equal(assignments.get('lens1.distribution'), '0.999f');
  assert.equal(assignments.get('lens1.speed'), '0.027205883f');
  assert.equal(assignments.get('lens1.drift'), '0.01632353f');
  assert.equal(assignments.get('source.lattice_cell_scale'), '2.420276f');
  assert.equal(assignments.get('modulation.modulation.sine.frequency'), '1.0f');
  assert.equal(assignments.has('value.edge_width'), false);
  source.descriptor.chain[0].label = 'whorl';
  source.descriptor.chain[0].resource_key = 'spin';
  const renamed = Object.fromEntries(Object.entries(values).map(([id, value]) =>
    [id.replace(/^lens1\./, 'whorl.'), value]));
  assert.equal(presetAssignments(source, renamed).get('spin.twist_angle'), '-10.259509f');
});


test('every catalog operator retains all scalar values and topology selections', () => {
  let operators = 0;
  let selectors = 0;
  for (const operator of catalog.operators) {
    const slot = { label: 'authored', operator: operator.id, resource_key: 'instance' };
    const authored = { descriptor: { chain: [slot] } };
    const values = Object.fromEntries(operator.params.map(field => [`authored.${field.id}`, field.default]));
    const assignments = presetAssignments(authored, values);
    for (const field of operator.params) {
      const id = `authored.${field.id}`;
      if (field.topology) {
        for (const choice of field.values) {
          const chosen = { ...values, [id]: choice };
          assert.equal(topologySelections(authored, chosen).get(id), choice, `${operator.id}/${id}/${choice}`);
          assert.doesNotThrow(() => presetAssignments(authored, chosen));
          ++selectors;
        }
      } else if (!(operator.id === 'sphere.rotate.v2' && field.id === 'spin-speed') && field.id !== 'edge-width') {
        assert.ok(assignments.size > 0, `${operator.id}/${id}`);
      }
    }
    ++operators;
  }
  assert.equal(operators, 44);
  assert.ok(selectors > 100);
});

test('operator bindings preserve authored resource keys and native member names', () => {
  const authored = { descriptor: { chain: [
    { label: 'fractal', operator: 'sample.fractal.v2' },
    { label: 'map', operator: 'project.s2.v1' },
    { label: 'flow', operator: 'sphere.displace.potential.v1', resource_key: 'velocity' },
    { label: 'slide', operator: 'warp.affine.v3', resource_key: 'plane' },
  ] } };
  assert.deepEqual([...presetAssignments(authored, {
    'fractal.fractal-scale': 2, 'fractal.fractal-iterations': 8,
    'fractal.julia-real': 0.5, 'fractal.julia-imaginary': -0.25,
    'fractal.fractal-contours': 3, 'map.coordinate-scale': 4,
    'flow.term-3-amplitude': 0.75, 'flow.term-3-phase': 0.125,
    'slide.lattice-period': 1.5,
  })], [
    ['source.scale', '2.0f'], ['source.iterations', '8.0f'],
    ['source.julia_re', '0.5f'], ['source.julia_im', '-0.25f'],
    ['source.contours', '3.0f'], ['projection.coordinate_scale', '4.0f'],
    ['velocity.terms[2].amplitude', '0.75f'], ['velocity.terms[2].phase', '0.125f'],
    ['plane.lattice_period', '1.5f'],
  ]);
});

test('authored topology emits each runtime catalog choice as its stored selector', () => {
  let checked = 0;
  for (const operator of catalog.operators) {
    if (operator.id === 'color.modulate.v1') continue;
    const authored = { descriptor: { chain: [{ label: 'node', operator: operator.id, resource_key: 'native' }] } };
    for (const field of operator.params.filter(field => field.topology)) {
      if (field.id === 'palette-mapping') continue;
      for (const [index, choice] of field.values.entries()) {
        const result = [...presetAssignments(authored, { [`node.${field.id}`]: choice }, { authoredTopology: true })];
        assert.equal(result.length, 1, `${operator.id}/${field.id}/${choice}`);
        if (field.id === 'brightness-envelope')
          assert.deepEqual(result[0], ['native.brightness_envelope', `Pullback::Color::BrightnessEnvelope::${choice.toUpperCase()}`]);
        else
          assert.equal(result[0][1], String(index), `${operator.id}/${field.id}/${choice}`);
        ++checked;
      }
    }
  }
  assert.ok(checked > 100);
});

test('camera spin binding follows its authored label', async () => {
  const authored = JSON.parse(await readFile(new URL('../patterns/grid_worms.shader.json', import.meta.url), 'utf8'));
  const rotate = authored.descriptor.chain.find(slot => slot.operator === 'sphere.rotate.v2');
  const original = rotate.label;
  rotate.label = 'view';
  for (const parameter of authored.descriptor.parameters)
    if (parameter.id.startsWith(`${original}.`)) parameter.id = parameter.id.replace(`${original}.`, 'view.');
  authored.descriptor.serialization.fields = authored.descriptor.serialization.fields.map(id => id.replace(`${original}.`, 'view.'));
  for (const preset of authored.preset_bank.presets) {
    for (const id of Object.keys(preset.values).filter(id => id.startsWith(`${original}.`))) {
      preset.values[id.replace(`${original}.`, 'view.')] = preset.values[id];
      delete preset.values[id];
    }
  }
  const assignments = presetAssignments(authored, { 'view.wander': 0.25 });
  assert.equal(assignments.get('projection.camera_wander'), '0.25f');
});

test('authored affine periods and sample coverage retain their independent storage', () => {
  const authored = { descriptor: { chain: [
    { label: 'warp', operator: 'warp.affine.v3' },
    { label: 'sample', operator: 'sample.lattice.v2' },
  ] } };
  const assignments = presetAssignments(authored, {
    'warp.lattice-period': 0.5, 'sample.lattice-cell-scale': 2,
    'sample.edge-width': 0.25, 'sample.weight-mode': 'projection',
    'sample.coverage-mode': 'edge-fade',
  }, { authoredTopology: true });
  assert.equal(assignments.get('outer_warp.lattice_period'), '0.5f');
  assert.equal(assignments.get('source.edge_width'), '0.25f');
  assert.equal(assignments.get('source.coverage_mode'), '3');
  assert.equal(assignments.get('source.weight_mode'), '1');
  assert.ok(!assignments.has('value.edge_width'));
});

test('authored rotations retain independent keys and per-preset spin speeds', async () => {
  const chain = { descriptor: { chain: [
    { label: 'near', operator: 'sphere.rotate.v2' },
    { label: 'far', operator: 'sphere.rotate.v2', resource_key: 'world' },
  ] } };
  assert.deepEqual([...presetAssignments(chain, {
    'near.wander': 0.25, 'near.spin-speed': 0.03125,
    'far.wander': 0.75, 'far.spin-speed': 0.015625,
  }, { authoredTopology: true })], [
    ['near.wander', '0.25f'], ['near.spin_speed', '0.03125f'],
    ['world.wander', '0.75f'], ['world.spin_speed', '0.015625f'],
  ]);
  const authored = JSON.parse(await readFile(new URL('../patterns/kaleidoscope_mandala.shader.json', import.meta.url), 'utf8'));
  const operator = catalog.operators.find(operator => operator.id === 'sphere.rotate.v2');
  authored.descriptor.parameters.push(declarationFromCatalogField('camera', operator.params.find(field => field.id === 'spin-speed'), operator.id));
  authored.descriptor.serialization.fields.push('camera.spin-speed');
  authored.preset_bank.presets[0].values['camera.spin-speed'] = 0;
  authored.preset_bank.presets[1].values['camera.spin-speed'] = 0.03125;
  const compiled = compileShaderDocument(authored, { catalog });
  assert.equal(compiled.status, 'VALID');
  const generated = generatedSections(compiled, undefined, { authoredTopology: true });
  assert.match(generated.params, /get<"camera">\(\)\.spin_speed = 0\.03125f;/);
  assert.doesNotMatch(generated.identity, /CAMERA_SPIN_RATE/);
  assert.throws(() => generatedSections(compiled), /camera spin must be uniform/);
});

test('periodic affine keeps explicit repeat units and authored winding controls', () => {
  const authored = { descriptor: { chain: [
    { label: 'repeat', resource_key: 'scroll', operator: 'warp.affine-periodic.v1' },
    { label: 'sample', operator: 'sample.lattice.v2' },
  ] } };
  for (const authoredTopology of [false, true]) {
    const assignments = presetAssignments(authored, {
      'repeat.lattice-period': 0.5, 'repeat.translation-x': -1.5,
      'repeat.translation-y': 2.25, 'sample.lattice-cell-scale': 2,
    }, { authoredTopology });
    assert.equal(assignments.get('scroll.lattice_period'), '0.5f');
    assert.equal(assignments.get('scroll.translation_x'), '-1.5f');
    assert.equal(assignments.get('scroll.translation_y'), '2.25f');
  }
});
