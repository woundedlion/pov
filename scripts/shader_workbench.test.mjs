import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile, readdir } from 'node:fs/promises';
import {
  ShaderDocumentError,
  canonicalPresetBank,
  classifyExport,
  compileShaderDocument,
  evaluateTransition,
  exportShaderDocumentJson,
  interpolateNormalizedGroup,
  interpolateValue,
  parseShaderDocument,
  stableStringify,
  validateShaderDocument,
} from './shader_workbench.mjs';
import { sha256Hex } from './sha256.mjs';

// scripts/shader_workbench.mjs, scripts/sha256.mjs and
// scripts/engine_catalog.json are verbatim mirrors of the daydream
// repository's shader/ copies (the engine install ships them there — see
// CMakeLists.txt). CI has no daydream checkout to diff against, so each mirror
// is pinned by the SHA-256 of its LF bytes instead: drift is a deliberate
// re-pin, updated together with the daydream master commit it mirrors.
// Mirrors daydream master e87bb510778c73597e74ffadf8a0bcc307a2396e.
// engine_catalog.json states the wasm32 operator ABI, the one the browser
// workbench's budget math models. tests/data/shader_chain_catalog.json is a
// separate catalog stating the native ABI the C++ suite pins. Their
// prepared-block sizes and alignments disagree by construction: the two files
// are not to be reconciled, and copying either over the other retargets a
// consumer's budget math.
const MIRROR_PINS = {
  'shader_workbench.mjs':
    'f31b0087ab00a3c34f578d5a77a490dec42a0f028275f6d01f772f84d3cb24c7',
  'sha256.mjs':
    '046a83178da02898524d3743ad3aa80ec91f719ea9c1b2e9f26912afba71015a',
  'engine_catalog.json':
    'd98634430ff8a0c6483130a8bcd1a9101b8ba92bd28469d3a69d4dcad14739ae',
};

const lf = (text) => text.replaceAll('\r\n', '\n');

const CATALOG = JSON.parse(
  await readFile(new URL('./engine_catalog.json', import.meta.url), 'utf8'));
const EXAMPLE = lf(await readFile(
  new URL('../patterns/example.shader.json', import.meta.url), 'utf8'));

/** @returns {Object} A fresh parse of the v2 example document, safe to mutate. */
const example = () => JSON.parse(EXAMPLE);

/** Compiles with the mirrored catalog. @param {*} source @param {Object} [options] */
const compile = (source, options = {}) =>
  compileShaderDocument(source, { catalog: CATALOG, ...options });

/** Validates with the mirrored catalog. @param {Object} document */
const validate = (document) =>
  validateShaderDocument(document, { catalog: CATALOG });

test('browser-compatible SHA-256 matches the published vectors', () => {
  assert.equal(sha256Hex(''),
    'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855');
  assert.equal(sha256Hex('abc'),
    'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad');
});

test('every mirrored compiler file matches its daydream pin', async () => {
  for (const [name, digest] of Object.entries(MIRROR_PINS)) {
    const text = await readFile(new URL(name, import.meta.url), 'utf8');
    assert.equal(sha256Hex(lf(text)), digest,
      `scripts/${name} drifted from the pinned daydream mirror`);
  }
});

// The native suite golden-pins tests/data/shader_chain_catalog.json from an
// LP64 host build; the mirror above is the same emitter's output from the
// wasm32 module, and is what an editor budgets arena bytes against. A
// pointer-bearing `prepared` block is wider under LP64, so the two disagree
// there by construction. This holds them to disagreeing about nothing else, so
// a golden regenerated on an unrelated host cannot re-pin quietly.
const POINTER_WIDENED_OPERATORS = 7;

test('the native golden and the wasm mirror differ only in pointer-block width', async () => {
  const golden = JSON.parse(await readFile(
    new URL('../tests/data/shader_chain_catalog.json', import.meta.url), 'utf8'));
  assert.equal(golden.catalog_version, CATALOG.catalog_version);
  assert.deepEqual(golden.budgets, CATALOG.budgets);
  assert.deepEqual(golden.carriers, CATALOG.carriers);
  assert.deepEqual(golden.operators.map((operator) => operator.id),
    CATALOG.operators.map((operator) => operator.id));
  let widened = 0;
  let widest = 0;
  golden.operators.forEach((native, index) => {
    const wasm = CATALOG.operators[index];
    assert.deepEqual({ ...wasm, blocks: null }, { ...native, blocks: null },
      `${native.id} diverges outside its block layout`);
    assert.deepEqual(wasm.blocks.param, native.blocks.param, `${native.id} param block`);
    assert.deepEqual(wasm.blocks.state, native.blocks.state, `${native.id} state block`);
    assert.ok(wasm.blocks.prepared.size <= native.blocks.prepared.size,
      `${native.id} prepared block is wider on wasm32 than on LP64`);
    if (wasm.blocks.prepared.size === native.blocks.prepared.size
      && wasm.blocks.prepared.align === native.blocks.prepared.align) return;
    widened += 1;
    widest = Math.max(widest, native.blocks.prepared.size - wasm.blocks.prepared.size);
    assert.equal(native.blocks.prepared.align, 8,
      `${native.id} prepared block diverges without carrying an LP64 pointer`);
    assert.equal(wasm.blocks.prepared.align, 4,
      `${native.id} prepared block diverges without narrowing to wasm32`);
  });
  assert.equal(widened, POINTER_WIDENED_OPERATORS,
    'the pointer-bearing operator set moved; re-pin it deliberately');
  // Worst case: a full-length chain of the widest divergent operator.
  assert.ok(widest * CATALOG.budgets.max_chain_ops < CATALOG.budgets.arena_bytes,
    'the two accountings can disagree by a whole arena');
});

test('every promoted shader document matches its compiled effect identity', async () => {
  const migration = JSON.parse(await readFile(
    new URL('../patterns/shaderball_migration.json', import.meta.url), 'utf8'));
  const headers = {
    'alien-ocean': 'AlienOcean.h',
    'contour-lattice': 'ContourLattice.h',
    'cosmic-eyeball': 'CosmicEyeball.h',
    'curl-facets': 'CurlFacets.h',
    'curl-lattice': 'CurlLattice.h',
    'equator-grid': 'EquatorGrid.h',
    'facet-grid': 'FacetGrid.h',
    'facet-wave': 'FacetWave.h',
    'glitch-grid': 'GlitchGrid.h',
    'hex-wave': 'HexWave.h',
    'kaleido-wave': 'KaleidoWave.h',
    'mobius-grid': 'MobiusGrid.h',
    'prism-lattice': 'PrismLattice.h',
    'prism-spiral': 'PrismSpiral.h',
    'signal-weave': 'SignalWeave.h',
    'vector-facets': 'VectorFacets.h',
  };
  assert.deepEqual(Object.keys(migration.source_documents).sort(),
    migration.product_group.children.map((child) => child.effect_id).sort());
  for (const [effectId, documentName] of Object.entries(migration.source_documents)) {
    const documentSource = await readFile(
      new URL(`../patterns/${documentName}`, import.meta.url), 'utf8');
    const header = await readFile(
      new URL(`../effects/${headers[effectId]}`, import.meta.url), 'utf8');
    const compiled = compile(parseShaderDocument(documentSource));
    assert.equal(compiled.status, 'VALID', effectId);
    assert.equal(compiled.document.effect_id, effectId);
    assert.ok(header.includes(`"${compiled.descriptor_digest}"`),
      `${effectId} descriptor digest is stale`);
    assert.ok(header.includes(`"${compiled.preset_bank_digest}"`),
      `${effectId} preset-bank digest is stale`);
    // The digests above cover the descriptor and preset values. Runtime timing
    // is fixed-effect policy, but PRESET_IDS must preserve the document's
    // generated order because preset_params() indexes it by position.
    const choreography = compiled.document.preset_bank.choreography;
    const declaredIds = header.match(
      /std::array<std::string_view,\s*(\d+)>\s+PRESET_IDS\{([^}]*)\}/);
    assert.ok(declaredIds, `${effectId} declares no PRESET_IDS array`);
    assert.deepEqual(
      [...declaredIds[2].matchAll(/"([^"]*)"/g)].map((match) => match[1]),
      choreography.generated_order,
      `${effectId} PRESET_IDS is not choreography.generated_order, in order`);
    assert.equal(Number(declaredIds[1]), choreography.generated_order.length,
      `${effectId} PRESET_IDS extent does not match the preset count`);
    const destinations = migration.destinations.filter((entry) => entry.effect_id === effectId);
    const presetIds = new Set(compiled.document.preset_bank.presets
      .map((preset) => preset.preset_id));
    for (const destination of destinations)
      assert.equal(presetIds.has(destination.preset_id), true,
        `${effectId}/${destination.preset_id} is missing from its source document`);
  }
});

test('every patterns/*.shader.json compiles and is accounted for', async () => {
  const directory = new URL('../patterns/', import.meta.url);
  const documents = (await readdir(directory))
    .filter((name) => name.endsWith('.shader.json')).sort();
  const migration = JSON.parse(await readFile(
    new URL('shaderball_migration.json', directory), 'utf8'));
  const promoted = new Set(Object.values(migration.source_documents));
  // The sample a contributor copies is not in the migration manifest, so the
  // promoted-document gate above never reaches it.
  assert.deepEqual(documents.filter((name) => !promoted.has(name)),
    ['example.shader.json']);
  for (const name of documents) {
    const compiled = compile(parseShaderDocument(
      await readFile(new URL(name, directory), 'utf8')));
    assert.equal(compiled.status, 'VALID', name);
  }
});

test('strict parsing rejects duplicate keys after NFC normalization', () => {
  assert.throws(
    () => parseShaderDocument('{"e\u0301":1,"é":2}'),
    (error) => error instanceof ShaderDocumentError && error.code === 'DUPLICATE_KEY',
  );
});

test('strict parsing enforces byte, depth, BOM, and finite-number bounds', () => {
  assert.throws(() => parseShaderDocument('\ufeff{}'), /byte-order mark/u);
  assert.throws(() => parseShaderDocument('{"value":1e999}'), /finite/u);
  assert.throws(() => parseShaderDocument('[[[]]]', { depth: 1 }), /nesting limit/u);
  assert.throws(() => parseShaderDocument('{"long":"abcd"}', { bytes: 4 }), /byte limit/u);
});

test('strict parsing skips only RFC 8259 whitespace', () => {
  for (const gap of ['\u00a0', '\u2028', '\u3000']) {
    assert.throws(() => parseShaderDocument(`${gap}{}`), ShaderDocumentError);
    assert.throws(() => parseShaderDocument(`{}${gap}`), ShaderDocumentError);
  }
});

test('the example chain document validates against the catalog', () => {
  assert.deepEqual(validate(example()), []);
});

test('unknown semantic fields are rejected', () => {
  const document = example();
  document.descriptor.chain[0].surprise = true;
  assert.throws(
    () => validate(document),
    (error) => error.code === 'UNKNOWN_FIELD' && error.path.endsWith('.surprise'),
  );
});

// A "__proto__" key has to land in the object: run through the prototype setter
// it leaves no unknown field to report and no bytes in the canonical form, so
// two different sources take one identity.
test('a __proto__ key is an ordinary field, not a prototype write', () => {
  const bare = parseShaderDocument('{"__proto__":{"schema_version":2}}');
  assert.deepEqual(Object.keys(bare), ['__proto__']);
  assert.equal('schema_version' in bare, false);

  const poisoned = compile(EXAMPLE.replace('"descriptor": {', '"descriptor": {"__proto__": {},'));
  assert.equal(poisoned.status, 'INVALID');
  assert.deepEqual(poisoned.diagnostics.map((entry) => [entry.code, entry.path]),
    [['UNKNOWN_FIELD', '$.descriptor.__proto__']]);

  // Metadata takes any key, so this document is accepted, and must stay distinct.
  const carried = compile(EXAMPLE.replace('"study_metadata": {', '"study_metadata": {"__proto__": 1,'));
  assert.equal(carried.status, 'VALID');
  assert.ok(exportShaderDocumentJson(carried.document).includes('"__proto__": 1'));
  assert.notEqual(exportShaderDocumentJson(carried.document),
    exportShaderDocumentJson(compile(EXAMPLE).document));
});

// The semantic phase collects and continues: the unknown operator and every
// parameter orphaned by it come back in one report.
test('an unknown operator and its orphaned parameters report together', () => {
  const document = example();
  document.descriptor.chain[2].operator = 'sample.future.v2';
  const diagnostics = validate(document);
  const codes = new Set(diagnostics.map((diagnostic) => diagnostic.code));
  assert.ok(codes.has('UNKNOWN_OPERATOR'));
  assert.ok(codes.has('UNBOUND_PARAMETER'));
  assert.equal(compile(document).status, 'INVALID');
});

test('chain carrier legality distinguishes order from mismatch', () => {
  const document = example();
  const chain = document.descriptor.chain;
  [chain[1], chain[2]] = [chain[2], chain[1]];
  const codes = new Set(validate(document).map((diagnostic) => diagnostic.code));
  assert.ok(codes.has('FAMILY_ORDER'));
  assert.ok(codes.has('CARRIER_MISMATCH'));
});

test('the descriptor digest survives reordering but not a label rename', () => {
  const baseline = compile(example());
  assert.equal(baseline.status, 'VALID');
  const shuffled = example();
  shuffled.descriptor.parameters.reverse();
  shuffled.preset_bank.presets.reverse();
  shuffled.preset_bank.edges.reverse();
  const compiled = compile(shuffled);
  assert.equal(compiled.status, 'VALID');
  assert.equal(compiled.descriptor_digest, baseline.descriptor_digest);
  assert.equal(compiled.preset_bank_digest, baseline.preset_bank_digest);

  // The chain labels and their order are digest-bearing: renaming an instance
  // no parameter binds is still a different descriptor.
  const renamed = example();
  renamed.descriptor.chain[0].label = 'rig';
  const recompiled = compile(renamed);
  assert.equal(recompiled.status, 'VALID');
  assert.notEqual(recompiled.descriptor_digest, baseline.descriptor_digest);
});

test('serialization fields name every parameter once and do not order the digest', () => {
  const baseline = compile(example());
  assert.equal(baseline.status, 'VALID');
  const reversed = example();
  reversed.descriptor.serialization.fields.reverse();
  assert.equal(compile(reversed).descriptor_digest, baseline.descriptor_digest,
    'field order is a spelling, not an identity');

  const codes = (mutate) => {
    const document = example();
    mutate(document.descriptor.serialization.fields);
    return validate(document).map((diagnostic) => diagnostic.code);
  };
  assert.deepEqual(codes((fields) => fields.splice(0)), ['INVALID_SERIALIZATION_FIELDS']);
  assert.deepEqual(codes((fields) => fields.pop()), ['INVALID_SERIALIZATION_FIELDS']);
  assert.deepEqual(codes((fields) => { fields[1] = fields[0]; }),
    ['INVALID_SERIALIZATION_FIELDS'], 'a duplicate hides a missing parameter');
  assert.deepEqual(codes((fields) => fields.push('sample.ghost-field')),
    ['INVALID_SERIALIZATION_FIELDS']);
});

test('preset values and document metadata do not enter semantic identity', () => {
  const first = example();
  const second = example();
  second.document_id = 'another-study';
  second.study_metadata.notes = 'different';
  second.preset_bank.presets[0].values['sample.pattern-freq'] = 2;
  second.preset_bank.presets[1].values['sample.coverage-mode'] = 'weight';
  assert.equal(compile(first).descriptor_digest, compile(second).descriptor_digest);
  assert.notEqual(compile(first).preset_bank_digest, compile(second).preset_bank_digest);
});

test('preset-bank identity ignores record declaration order but preserves generated order', () => {
  const first = example();
  const second = example();
  second.preset_bank.presets.reverse();
  second.preset_bank.edges.reverse();
  assert.equal(stableStringify(canonicalPresetBank(first)),
    stableStringify(canonicalPresetBank(second)));
  second.preset_bank.choreography.generated_order.reverse();
  assert.notEqual(stableStringify(canonicalPresetBank(first)),
    stableStringify(canonicalPresetBank(second)));
});

// The daydream v1 example fixture, inlined: expanding it must reproduce the
// committed v2 example byte for byte, pinning expandV1Document as the single
// code path both schema generations share.
const V1_EXAMPLE = {
  schema_version: 1,
  catalog_version: 1,
  document_id: 'example-study',
  effect_id: null,
  descriptor: {
    graph: {
      nodes: [
        { label: 'camera', role: 'outer_camera', operator: 'pullback.outer_camera.v1' },
        { label: 'surface', role: 'surface_project', operator: 'pullback.surface_project.v1',
          policy: { lens: 'identity', projection: 'equirectangular' } },
        { label: 'warp', role: 'planar_warp', operator: 'pullback.planar_warp.v1',
          policy: { outer: 'identity' } },
        { label: 'pattern', role: 'source', operator: 'pullback.source.v1',
          policy: { source: 'grid' } },
        { label: 'transfer', role: 'material', operator: 'pullback.material.v1',
          policy: { weight: 'projection', transfer: 'linear', coverage: 'opaque' } },
        { label: 'palette', role: 'color', operator: 'pullback.color.v1',
          policy: { color: 'generated_palette' } },
      ],
      edges: [
        { from: 'camera', to: 'surface' },
        { from: 'surface', to: 'warp' },
        { from: 'warp', to: 'pattern' },
        { from: 'pattern', to: 'transfer' },
        { from: 'transfer', to: 'palette' },
      ],
    },
    parameters: [
      {
        id: 'pattern-freq', binding: 'source.pattern-freq', classification: 'preset',
        storage: 'binary32', unit: 'ratio', domain: { minimum: 0.01, maximum: 10 },
        interpolation: { kind: 'LOG_POSITIVE' }, default: 1,
      },
      {
        id: 'central-meridian', binding: 'projection.central-meridian', classification: 'preset',
        storage: 'binary32', unit: 'radian',
        domain: { minimum: 0, maximum: 6.2831854820251465 },
        interpolation: { kind: 'SHORTEST_PERIODIC', period: 6.2831854820251465 },
        default: 0,
      },
    ],
    path_policies: [{ id: 'parallel', kind: 'PARALLEL' }],
    clocks: [{ id: 'source-clock', kind: 'frame-clock', settings: { wrap: 1 } }],
    preparation: [{ id: 'frame', kind: 'prepare-frame' }],
    resources: [],
    serialization: { schema_version: 1, fields: ['pattern-freq', 'central-meridian'] },
    approximation: [],
    handoff: { policy: 'reset' },
  },
  preset_bank: {
    schema_version: 1,
    presets: [
      { preset_id: 'calm', display_name: 'Calm',
        values: { 'pattern-freq': 1, 'central-meridian': 6 } },
      { preset_id: 'fast', display_name: 'Fast',
        values: { 'pattern-freq': 4, 'central-meridian': 0.2 } },
    ],
    edges: [
      { from: 'calm', to: 'fast', path_policy: 'parallel', easing: 'EASE_IN_OUT_SIN', duration: 120 },
      { from: 'fast', to: 'calm', path_policy: 'parallel', easing: 'EASE_IN_OUT_SIN', duration: 120 },
    ],
    absent_edge_fallback: {
      manual: 'SNAP', automatic: 'REJECT', synchronized: 'REJECT',
      restore: 'SNAP', authoring: 'SNAP',
    },
    choreography: { generated_order: ['calm', 'fast'], dwell: { calm: 600, fast: 600 } },
  },
  study_metadata: { notes: 'Example authoring document' },
};

test('a v1 document expands to the committed v2 example byte for byte', () => {
  const compiled = compile(structuredClone(V1_EXAMPLE));
  assert.equal(compiled.status, 'VALID');
  assert.equal(compiled.parameter_ids['pattern-freq'], 'sample.pattern-freq');
  assert.equal(compiled.parameter_ids['central-meridian'], 'project.central-meridian');
  assert.equal(exportShaderDocumentJson(compiled.document), EXAMPLE);
  assert.equal(compiled.descriptor_digest, compile(example()).descriptor_digest);
});

test('export classification compares exact descriptors after the digest', () => {
  const compiled = compile(example());
  const registry = { effects: [{
    effect_id: 'curl-lattice',
    descriptor_digest: compiled.descriptor_digest,
    descriptor: compiled.descriptor,
    capability_profiles: ['wasm-authoring'],
  }] };
  assert.deepEqual(classifyExport(compiled, registry, 'wasm-authoring'),
    { kind: 'ADD_PRESET_CANDIDATE', effect_id: 'curl-lattice' });
  registry.effects[0].descriptor = { ...compiled.descriptor, serialization: { schema_version: 2, fields: [] } };
  assert.equal(classifyExport(compiled, registry, 'wasm-authoring').kind, 'CREATE_EFFECT_CANDIDATE');
});

test('known but target-unavailable effects remain distinguishable', () => {
  const compiled = compile(example());
  const result = classifyExport(compiled, { effects: [{
    effect_id: 'curl-lattice', descriptor_digest: compiled.descriptor_digest,
    descriptor: compiled.descriptor, capability_profiles: ['wasm-authoring'],
  }] }, 'teensy-shipping');
  assert.equal(result.kind, 'REJECTED');
  assert.equal(result.effect_id, 'curl-lattice');
  assert.equal(result.diagnostics[0].code, 'KNOWN_UNAVAILABLE');
});

test('linear and log interpolation preserve exact stored endpoints', () => {
  const scale = example().descriptor.parameters
    .find((parameter) => parameter.id === 'sample.pattern-freq');
  const linear = structuredClone(scale);
  linear.interpolation = { kind: 'LINEAR' };
  assert.equal(interpolateValue(linear, 1, 4, -1), Math.fround(1));
  assert.equal(interpolateValue(linear, 1, 4, 2), Math.fround(4));
  assert.equal(interpolateValue(scale, 1, 4, 0.5), Math.fround(2));
});

test('mixed enum parameters expose blend state between distinct endpoints', () => {
  const mapping = example().descriptor.parameters
    .find((parameter) => parameter.id === 'sample.coverage-mode');
  assert.deepEqual(interpolateValue(mapping, 'none', 'weight', 0.25),
    { from: 'none', to: 'weight', mix: Math.fround(0.25) });
  assert.equal(interpolateValue(mapping, 'none', 'weight', 0), 'none');
  assert.equal(interpolateValue(mapping, 'none', 'weight', 1), 'weight');
  assert.equal(interpolateValue(mapping, 'edge-fade', 'edge-fade', 0.5), 'edge-fade');
});

test('periodic interpolation uses the negative half-period tie', () => {
  const phase = {
    id: 'phase', classification: 'preset', storage: 'binary32', unit: 'turn',
    domain: { minimum: 0, maximum: 1 },
    interpolation: { kind: 'SHORTEST_PERIODIC', period: 1 }, default: 0,
  };
  assert.equal(interpolateValue(phase, 0, 0.5, 0.5), Math.fround(0.75));
  assert.ok(Math.abs(interpolateValue(phase, 0.9, 0.1, 0.5)) < 1e-6);
});

test('a transition evaluates its exact endpoint values', () => {
  const document = example();
  const source = evaluateTransition(document.descriptor, document.preset_bank, 'fast', 'calm', 0);
  const destination = evaluateTransition(document.descriptor, document.preset_bank, 'fast', 'calm', 120);
  assert.deepEqual(source.values, {
    'project.central-meridian': Math.fround(0.2),
    'sample.coverage-mode': 'none',
    'sample.pattern-freq': Math.fround(4),
    'sample.weight-mode': 'projection',
  });
  assert.deepEqual(destination.values, {
    'project.central-meridian': Math.fround(6),
    'sample.coverage-mode': 'none',
    'sample.pattern-freq': Math.fround(1),
    'sample.weight-mode': 'projection',
  });
});

test('staggered paths apply easing before ordered group scheduling', () => {
  const document = example();
  document.descriptor.path_policies[0] = {
    id: 'parallel', kind: 'STAGGERED_ORDERED',
    groups: ['sample.pattern-freq', 'project.central-meridian',
      'sample.weight-mode', 'sample.coverage-mode'],
  };
  const result = evaluateTransition(document.descriptor, document.preset_bank, 'calm', 'fast', 60);
  assert.equal(result.eased_progress, Math.fround(0.5));
  assert.equal(result.values['sample.pattern-freq'], Math.fround(4));
  assert.equal(result.values['project.central-meridian'], Math.fround(0.2));
  assert.equal(result.values['sample.coverage-mode'], 'none');
});

test('normalized interpolation evaluates complete groups and rejects antipodes', () => {
  const fields = ['x', 'y'].map((field) => ({
    id: field, domain: { minimum: -1, maximum: 1 },
    interpolation: { kind: 'NORMALIZED_LINEAR', group: 'axis' },
  }));
  assert.deepEqual(interpolateNormalizedGroup(fields, { x: 1, y: 0 }, { x: 0, y: 1 }, 0.5),
    { x: Math.fround(Math.SQRT1_2), y: Math.fround(Math.SQRT1_2) });
  assert.throws(
    () => interpolateNormalizedGroup(fields, { x: 1, y: 0 }, { x: -1, y: 0 }, 0.5),
    (error) => error.code === 'DEGENERATE_NORMALIZED_PATH',
  );
});

test('every promoted ShaderWorkbench preset has one stable migration destination', async () => {
  const migration = JSON.parse(await readFile(
    new URL('../patterns/shaderball_migration.json', import.meta.url), 'utf8'));
  assert.equal(migration.legacy_alias, 'ShaderBall');
  assert.equal(migration.authoring_effect, 'Shader');
  assert.deepEqual(migration.retired_legacy_presets, [4]);
  assert.deepEqual(migration.destinations.map((entry) => entry.legacy_preset),
    Array.from({ length: 24 }, (_, index) => index).filter((index) => index !== 4));
  assert.equal(new Set(migration.destinations
    .map((entry) => `${entry.effect_id}/${entry.preset_id}`)).size, 23);
  assert.equal(migration.product_group.children
    .reduce((total, child) => total + child.seconds, 0), 120);
  const childIds = new Set(migration.product_group.children
    .map((child) => child.effect_id));
  for (const destination of migration.destinations) {
    assert.equal(childIds.has(destination.effect_id), true,
      `${destination.effect_id} is missing from product discovery`);
  }
});
