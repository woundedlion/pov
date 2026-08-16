import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import {
  ShaderDocumentError,
  canonicalDescriptor,
  canonicalPresetBank,
  classifyExport,
  compileShaderDocument,
  evaluateTransition,
  interpolateNormalizedGroup,
  interpolateValue,
  parseShaderDocument,
  stableStringify,
  validateShaderDocument,
} from './shader_workbench.mjs';

const makeDocument = () => ({
  schema_version: 1,
  catalog_version: 1,
  document_id: 'curl-study',
  effect_id: null,
  descriptor: {
    graph: {
      nodes: [
        { label: 'camera', role: 'outer_camera', operator: 'pullback.outer_camera.v1' },
        { label: 'surface', role: 'surface_project', operator: 'pullback.surface_project.v1', policy: { lens: 'identity' } },
        { label: 'warp', role: 'planar_warp', operator: 'pullback.planar_warp.v1' },
        { label: 'pattern', role: 'source', operator: 'pullback.source.v1' },
        { label: 'transfer', role: 'material', operator: 'pullback.material.v1' },
        { label: 'palette', role: 'color', operator: 'pullback.color.v1' },
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
        id: 'scale', binding: 'source.scale', classification: 'preset', storage: 'binary32',
        unit: 'ratio', domain: { minimum: 0.01, maximum: 10 },
        interpolation: { kind: 'LOG_POSITIVE' }, default: 1,
      },
      {
        id: 'phase', binding: 'source.phase', classification: 'preset', storage: 'binary32',
        unit: 'turn', domain: { minimum: 0, maximum: 1 },
        interpolation: { kind: 'SHORTEST_PERIODIC', period: 1 }, default: 0,
      },
    ],
    path_policies: [{ id: 'parallel', kind: 'PARALLEL' }],
    clocks: [{ id: 'source-clock', kind: 'frame-clock', settings: { wrap: 1 } }],
    preparation: [{ id: 'frame', kind: 'prepare-frame' }],
    resources: [],
    serialization: { schema_version: 1, fields: ['scale', 'phase'] },
    approximation: [],
    handoff: { policy: 'reset' },
  },
  preset_bank: {
    schema_version: 1,
    presets: [
      { preset_id: 'calm', display_name: 'Calm', values: { scale: 1, phase: 0.9 } },
      { preset_id: 'fast', display_name: 'Fast', values: { scale: 4, phase: 0.1 } },
    ],
    edges: [
      { from: 'calm', to: 'fast', path_policy: 'parallel', easing: 'LINEAR', duration: 4 },
      { from: 'fast', to: 'calm', path_policy: 'parallel', easing: 'EASE_IN_OUT_SIN', duration: 1 },
    ],
    absent_edge_fallback: {
      manual: 'SNAP', automatic: 'REJECT', synchronized: 'REJECT', restore: 'SNAP', authoring: 'SNAP',
    },
    choreography: { generated_order: ['calm', 'fast'], dwell: { calm: 120, fast: 120 } },
  },
  study_metadata: { notes: 'authoring only' },
});

test('shader documents match their compiled effect identities', async () => {
  for (const [documentName, headerName] of [
    ['curl_lattice.shader.json', 'CurlLattice.h'],
    ['facet_grid.shader.json', 'FacetGrid.h'],
  ]) {
    const documentSource = await readFile(
      new URL(`../patterns/${documentName}`, import.meta.url), 'utf8');
    const header = await readFile(
      new URL(`../effects/${headerName}`, import.meta.url), 'utf8');
    const compiled = compileShaderDocument(parseShaderDocument(documentSource));
    const descriptor = header.match(/DESCRIPTOR_DIGEST\s*=\s*\n\s*"([0-9a-f]{64})"/u);
    const bank = header.match(/PRESET_BANK_DIGEST\s*=\s*\n\s*"([0-9a-f]{64})"/u);
    assert.ok(descriptor);
    assert.ok(bank);
    assert.equal(descriptor[1], compiled.descriptor_digest);
    assert.equal(bank[1], compiled.preset_bank_digest);
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

test('the linear document validates with exhaustive parameter values', () => {
  assert.deepEqual(validateShaderDocument(makeDocument()), []);
});

test('unknown semantic fields are rejected', () => {
  const document = makeDocument();
  document.descriptor.graph.nodes[0].surprise = true;
  assert.throws(
    () => validateShaderDocument(document),
    (error) => error.code === 'UNKNOWN_FIELD' && error.path.endsWith('.surprise'),
  );
});

test('invalid graph topology is rejected independently of declaration order', () => {
  const document = makeDocument();
  document.descriptor.graph.edges[2] = { from: 'warp', to: 'transfer' };
  assert.throws(
    () => validateShaderDocument(document),
    (error) => error.code === 'INVALID_LINEAR_GRAPH',
  );
});

test('canonical descriptors ignore labels and declaration order', () => {
  const first = makeDocument();
  const second = structuredClone(first);
  second.descriptor.graph.nodes.reverse();
  second.descriptor.graph.edges.reverse();
  for (const node of second.descriptor.graph.nodes) node.label = `renamed-${node.role}`;
  second.descriptor.graph.edges = [
    { from: 'renamed-outer_camera', to: 'renamed-surface_project' },
    { from: 'renamed-surface_project', to: 'renamed-planar_warp' },
    { from: 'renamed-planar_warp', to: 'renamed-source' },
    { from: 'renamed-source', to: 'renamed-material' },
    { from: 'renamed-material', to: 'renamed-color' },
  ].reverse();
  second.descriptor.parameters.reverse();
  assert.equal(stableStringify(canonicalDescriptor(first)), stableStringify(canonicalDescriptor(second)));
  assert.equal(compileShaderDocument(first).descriptor_digest, compileShaderDocument(second).descriptor_digest);
});

test('preset values and document metadata do not enter semantic identity', () => {
  const first = makeDocument();
  const second = structuredClone(first);
  second.document_id = 'another-study';
  second.study_metadata.notes = 'different';
  second.preset_bank.presets[0].values.scale = 2;
  assert.equal(compileShaderDocument(first).descriptor_digest, compileShaderDocument(second).descriptor_digest);
  assert.notEqual(compileShaderDocument(first).preset_bank_digest, compileShaderDocument(second).preset_bank_digest);
});

test('preset-bank identity ignores record declaration order but preserves generated order', () => {
  const first = makeDocument();
  const second = structuredClone(first);
  second.preset_bank.presets.reverse();
  second.preset_bank.edges.reverse();
  assert.equal(stableStringify(canonicalPresetBank(first)), stableStringify(canonicalPresetBank(second)));
  second.preset_bank.choreography.generated_order.reverse();
  assert.notEqual(stableStringify(canonicalPresetBank(first)), stableStringify(canonicalPresetBank(second)));
});

test('an unavailable catalog operator is retained but blocks export', () => {
  const document = makeDocument();
  document.descriptor.graph.nodes[1].operator = 'pullback.future_branch.v1';
  const compiled = compileShaderDocument(document);
  assert.equal(compiled.status, 'VALID_BUT_UNSUPPORTED');
  assert.equal(compiled.diagnostics[0].code, 'OPERATOR_UNAVAILABLE');
});

test('export classification compares exact descriptors after the digest', () => {
  const compiled = compileShaderDocument(makeDocument());
  const registry = { effects: [{
    effect_id: 'curl-lattice',
    descriptor_digest: compiled.descriptor_digest,
    descriptor: compiled.descriptor,
    capability_profiles: ['wasm-authoring'],
  }] };
  assert.deepEqual(classifyExport(compiled, registry, 'wasm-authoring'),
    { kind: 'ADD_PRESET_CANDIDATE', effect_id: 'curl-lattice' });
  registry.effects[0].descriptor = { ...compiled.descriptor, handoff: { policy: 'exact' } };
  assert.equal(classifyExport(compiled, registry, 'wasm-authoring').kind, 'CREATE_EFFECT_CANDIDATE');
});

test('known but target-unavailable effects remain distinguishable', () => {
  const compiled = compileShaderDocument(makeDocument());
  const result = classifyExport(compiled, { effects: [{
    effect_id: 'curl-lattice', descriptor_digest: compiled.descriptor_digest,
    descriptor: compiled.descriptor, capability_profiles: ['wasm-authoring'],
  }] }, 'teensy-shipping');
  assert.equal(result.kind, 'REJECTED');
  assert.equal(result.effect_id, 'curl-lattice');
  assert.equal(result.diagnostics[0].code, 'KNOWN_UNAVAILABLE');
});

test('linear and log interpolation preserve exact stored endpoints', () => {
  const [scale] = makeDocument().descriptor.parameters;
  const linear = structuredClone(scale);
  linear.interpolation = { kind: 'LINEAR' };
  assert.equal(interpolateValue(linear, 1, 4, -1), Math.fround(1));
  assert.equal(interpolateValue(linear, 1, 4, 2), Math.fround(4));
  assert.equal(interpolateValue(scale, 1, 4, 0.5), Math.fround(2));
});

test('mixed enum parameters stay in one descriptor and expose blend state', () => {
  const first = makeDocument();
  const mapping = {
    id: 'palette-mapping', binding: 'color.palette-mapping',
    classification: 'preset', storage: 'enum8', unit: 'mapping',
    domain: { values: ['cup', 'bell', 'linear', 'reverse'] },
    interpolation: { kind: 'MIXED_ENUM' }, default: 'cup',
  };
  first.descriptor.parameters.push(mapping);
  first.descriptor.serialization.fields.push(mapping.id);
  first.preset_bank.presets[0].values[mapping.id] = 'cup';
  first.preset_bank.presets[1].values[mapping.id] = 'bell';
  const second = structuredClone(first);
  second.preset_bank.presets[1].values[mapping.id] = 'reverse';

  assert.equal(compileShaderDocument(first).descriptor_digest,
    compileShaderDocument(second).descriptor_digest);
  assert.notEqual(compileShaderDocument(first).preset_bank_digest,
    compileShaderDocument(second).preset_bank_digest);
  assert.deepEqual(interpolateValue(mapping, 'cup', 'bell', 0.25),
    { from: 'cup', to: 'bell', mix: Math.fround(0.25) });
  assert.equal(interpolateValue(mapping, 'cup', 'bell', 0), 'cup');
  assert.equal(interpolateValue(mapping, 'cup', 'bell', 1), 'bell');
});

test('periodic interpolation uses the negative half-period tie', () => {
  const phase = makeDocument().descriptor.parameters[1];
  assert.equal(interpolateValue(phase, 0, 0.5, 0.5), Math.fround(0.75));
  assert.ok(Math.abs(interpolateValue(phase, 0.9, 0.1, 0.5)) < 1e-6);
});

test('duration one evaluates source and destination on successive frames', () => {
  const document = makeDocument();
  const source = evaluateTransition(document.descriptor, document.preset_bank, 'fast', 'calm', 0);
  const destination = evaluateTransition(document.descriptor, document.preset_bank, 'fast', 'calm', 1);
  assert.deepEqual(source.values, { scale: Math.fround(4), phase: Math.fround(0.1) });
  assert.deepEqual(destination.values, { scale: Math.fround(1), phase: Math.fround(0.9) });
});

test('staggered paths apply easing before ordered group scheduling', () => {
  const document = makeDocument();
  document.descriptor.path_policies[0] = {
    id: 'parallel', kind: 'STAGGERED_ORDERED', groups: ['scale', 'phase'],
  };
  const result = evaluateTransition(document.descriptor, document.preset_bank, 'calm', 'fast', 1);
  assert.equal(result.values.scale, Math.fround(2));
  assert.equal(result.values.phase, Math.fround(0.9));
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

test('every ShaderBall preset has one stable migration destination', async () => {
  const migration = JSON.parse(await readFile(
    new URL('../patterns/shaderball_migration.json', import.meta.url), 'utf8'));
  assert.equal(migration.legacy_alias, 'ShaderBall');
  assert.equal(migration.authoring_effect, 'Shader');
  assert.deepEqual(migration.destinations.map((entry) => entry.legacy_preset),
    Array.from({ length: 19 }, (_, index) => index));
  assert.equal(new Set(migration.destinations
    .map((entry) => `${entry.effect_id}/${entry.preset_id}`)).size, 19);
  assert.equal(migration.product_group.children
    .reduce((total, child) => total + child.seconds, 0), 120);
  const childIds = new Set(migration.product_group.children
    .map((child) => child.effect_id));
  for (const destination of migration.destinations) {
    assert.equal(childIds.has(destination.effect_id), true,
      `${destination.effect_id} is missing from product discovery`);
  }
});
