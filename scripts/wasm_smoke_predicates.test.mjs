// Unit tests for the WASM smoke test's decisions (wasm_smoke_predicates.mjs).
//
// Run: node --test scripts/wasm_smoke_predicates.test.mjs

import test from 'node:test';
import assert from 'node:assert/strict';

import {
  BAKED_CONSTANT_IDS,
  DARK_BAND,
  DARK_EXEMPT,
  STACK_MAX_FILL,
  bakedTopologyFields,
  darknessExpectation,
  darknessProblems,
  engineControlNames,
  paramStreamProblems,
  promotedBindingProblems,
  stackCreepBudget,
} from './wasm_smoke_predicates.mjs';

const RES = [[288, 144], [96, 20]];
const KEYS = RES.map(([w, h]) => `${DARK_EXEMPT}@${w}x${h}`);

const sweep = (extra = {}) => ({
  frames: 3,
  resolutions: RES,
  sweptEffects: new Set([DARK_EXEMPT, 'Voronoi']),
  darkKeys: new Set(KEYS),
  ...extra,
});

test('a short window expects the exempt effect dark', () => {
  assert.equal(darknessExpectation(1), 'exempt-dark');
  assert.equal(darknessExpectation(DARK_BAND[0] - 1), 'exempt-dark');
});

test('the ambiguous band accepts either verdict', () => {
  assert.equal(darknessExpectation(DARK_BAND[0]), 'either');
  assert.equal(darknessExpectation(DARK_BAND[1] - 1), 'either');
});

test('a long window expects every effect lit', () => {
  assert.equal(darknessExpectation(DARK_BAND[1]), 'lit');
  assert.equal(darknessExpectation(1000), 'lit');
});

test('a short window accepts exactly the exempt effect being dark', () => {
  assert.deepEqual(darknessProblems(sweep()), []);
});

test('a short window fails when the exempt effect lit anyway', () => {
  const problems = darknessProblems(sweep({ darkKeys: new Set([KEYS[0]]) }));
  assert.equal(problems.length, 1);
  assert.match(problems[0], /lit after 3 frame\(s\)/);
});

test('any other dark effect fails at every window length', () => {
  for (const frames of [3, DARK_BAND[0], DARK_BAND[1]]) {
    const problems = darknessProblems(sweep({
      frames,
      darkKeys: new Set([`Voronoi@288x144`]),
    }));
    assert.ok(problems.some((p) => p.startsWith('Voronoi@288x144:')), `frames=${frames}`);
  }
});

test('the ambiguous band accepts the exempt effect dark or lit', () => {
  const frames = DARK_BAND[0];
  assert.deepEqual(darknessProblems(sweep({ frames })), []);
  assert.deepEqual(darknessProblems(sweep({ frames, darkKeys: new Set() })), []);
});

test('a long window rejects the exempt effect staying dark', () => {
  const problems = darknessProblems(sweep({ frames: DARK_BAND[1] }));
  assert.equal(problems.length, RES.length);
  assert.ok(problems.every((p) => /every pixel is zero/.test(p)));
});

test('an exemption naming no rendered effect is stale', () => {
  const problems = darknessProblems(sweep({
    frames: DARK_BAND[1],
    sweptEffects: new Set(['Voronoi']),
    darkKeys: new Set(),
  }));
  assert.equal(problems.length, 1);
  assert.match(problems[0], /the exemption is stale/);
});

test('the stack budget is the lower of the ceiling and the capacity fraction', () => {
  assert.equal(stackCreepBudget({ capacity: 1000 }, 2048), 1000 * STACK_MAX_FILL);
  assert.equal(stackCreepBudget({ capacity: 65536 }, 2048), 2048);
});

test('a missing or degenerate stack falls back to the ceiling', () => {
  assert.equal(stackCreepBudget(undefined, 2048), 2048);
  assert.equal(stackCreepBudget({ capacity: 0 }, 2048), 2048);
  assert.equal(stackCreepBudget({ capacity: -8 }, 2048), 2048);
});

const floatDef = (name, value, min = 0, max = 1) => ({ name, value, min, max });

test('matching parameter streams zip clean', () => {
  const defs = [floatDef('Speed', 0.25), floatDef('Scale', 3, 1, 8),
    { name: 'Wrap', value: true }];
  assert.deepEqual(paramStreamProblems(defs, [0.25, 3, 1]), []);
});

test('a non-array definition stream is the only problem reported', () => {
  assert.deepEqual(paramStreamProblems(undefined, []),
    ['getParameterDefinitions() did not return an array']);
});

test('a non-indexable value stream is the only problem reported', () => {
  const expected = ['getParamValues() did not return an indexable value stream'];
  assert.deepEqual(paramStreamProblems([floatDef('Speed', 0.25)], undefined),
    expected);
  assert.deepEqual(paramStreamProblems([floatDef('Speed', 0.25)], {}), expected);
});

test('a definition entry that is not an object is reported', () => {
  const problems = paramStreamProblems([null], [0.25]);
  assert.deepEqual(problems, ['param 0 is not a definition object']);
});

test('a length split is reported', () => {
  const problems = paramStreamProblems([floatDef('Speed', 0.25)], []);
  assert.ok(problems.some((p) => /param order seam drifted/.test(p)));
});

test('a transposition of two same-range floats is caught', () => {
  const defs = [floatDef('Speed', 0.25), floatDef('Scale', 0.75)];
  const problems = paramStreamProblems(defs, [0.75, 0.25]);
  assert.equal(problems.length, 2);
  assert.ok(problems.every((p) => /param order seam transposed/.test(p)));
});

test('a bool def is compared against the raw float it was collapsed from', () => {
  const defs = [{ name: 'Wrap', value: true }];
  assert.deepEqual(paramStreamProblems(defs, [0.9]), []);
  assert.deepEqual(paramStreamProblems(defs, [1]), []);
  const problems = paramStreamProblems(defs, [0.1]);
  assert.equal(problems.length, 1);
  assert.match(problems[0], /def bool true != value-stream 0.1 > 0.5/);
});

test('a bool def is exempt from the range check it carries no bounds for', () => {
  assert.deepEqual(paramStreamProblems([{ name: 'Wrap', value: false }], [0]), []);
});

test('a nameless parameter is reported', () => {
  const problems = paramStreamProblems([{ name: '', value: 0.5, min: 0, max: 1 }], [0.5]);
  assert.ok(problems.some((p) => p === 'param 0 has no name'));
});

test('a non-finite value-stream entry is reported', () => {
  const problems = paramStreamProblems([floatDef('Speed', 0.5)], [Number.NaN]);
  assert.ok(problems.some((p) => /is not finite/.test(p)));
});

test('an inverted or non-finite range is reported', () => {
  assert.ok(paramStreamProblems([floatDef('Speed', 0.5, 1, 0)], [0.5])
    .some((p) => /non-finite\/inverted range/.test(p)));
  assert.ok(paramStreamProblems([floatDef('Speed', 0.5, 0, Number.POSITIVE_INFINITY)], [0.5])
    .some((p) => /non-finite\/inverted range/.test(p)));
});

test('a value outside its range is reported, with a tolerance for rounding', () => {
  assert.ok(paramStreamProblems([floatDef('Speed', 2)], [2])
    .some((p) => /outside \[0, 1\]/.test(p)));
  // Within the relative epsilon the gate allows for a float round trip.
  assert.deepEqual(paramStreamProblems([floatDef('Speed', 1 + 1e-5)], [1 + 1e-5]), []);
});

const CONSTANT_ID = [...BAKED_CONSTANT_IDS][0];

const bindings = (extra = {}) => ({
  documents: [{
    document: 'patterns/ash_cloud.shader.json',
    effect: 'ash-cloud',
    parameterIds: ['camera.wander', 'surface.scale', 'lens.symmetry', CONSTANT_ID],
  }],
  controls: new Map([['ash-cloud', new Set(['Camera Wander', 'Surface Noise Scale'])]]),
  bakedFields: new Set(['symmetry']),
  ...extra,
});

test('a topology field and a baked constant need no control', () => {
  assert.deepEqual(promotedBindingProblems(bindings()), []);
});

test('an id that names no control is reported with what it tried', () => {
  const problems = promotedBindingProblems(bindings({
    documents: [{
      document: 'patterns/ash_cloud.shader.json',
      effect: 'ash-cloud',
      parameterIds: [CONSTANT_ID, 'sample.edge-width'],
    }],
  }));
  assert.equal(problems.length, 1);
  assert.match(problems[0], /"sample.edge-width" names no control on "ash-cloud"/);
  assert.match(problems[0], /tried "Edge Width"/);
});

test('a document naming an absent effect is one problem, not one per id', () => {
  const problems = promotedBindingProblems(bindings({ controls: new Map() }));
  assert.equal(problems.length, 2);
  assert.match(problems[0], /carries no effect "ash-cloud"/);
  assert.match(problems[1], /exemption is stale/);
});

test('a baked-constant exemption no document carries is stale', () => {
  const problems = promotedBindingProblems(bindings({
    documents: [{
      document: 'patterns/lattice_melt.shader.json',
      effect: 'ash-cloud',
      parameterIds: ['camera.wander'],
    }],
  }));
  assert.equal(problems.length, 1);
  assert.match(problems[0], /exemption is stale/);
});

test('the alias table resolves the label-derived engine names', () => {
  assert.deepEqual(engineControlNames('camera.wander'), ['Camera Wander']);
  assert.deepEqual(engineControlNames('surface.speed'), ['Surface Noise Speed']);
  assert.deepEqual(engineControlNames('sample.angle-speed'), ['Source Angle Speed']);
  assert.deepEqual(engineControlNames('cutout.cutout-threshold'), ['Cutout Threshold']);
  assert.deepEqual(engineControlNames('project.singularity-fade'), ['Singularity Fade']);
  assert.deepEqual(engineControlNames('palette-chroma'), ['Palette Chroma']);
});

test('a warp slot falls back from its slot name to its family name', () => {
  assert.deepEqual(engineControlNames('warp1.speed'), ['Planar Warp 1 Speed', 'Speed']);
  assert.deepEqual(engineControlNames('warp2.shear'), ['Planar Warp 2 Shear', 'Affine Shear']);
  assert.deepEqual(engineControlNames('warp2.cell-x'), ['Planar Warp 2 Cell X', 'Mirror Cell X']);
  assert.deepEqual(engineControlNames('warp1.radial-phase'),
    ['Planar Warp 1 Radial Phase', 'Polar Radial Phase']);
  assert.deepEqual(engineControlNames('warp1.field-angle'),
    ['Planar Warp 1 Field Angle', 'Warp Field Angle']);
});

test('the baked topology fields are the catalog flag less the live mapping', () => {
  const fields = bakedTopologyFields({
    operators: [{ params: [
      { id: 'symmetry', topology: true },
      { id: 'palette-mapping', topology: true },
      { id: 'lattice-radius' },
    ] }],
  });
  assert.deepEqual([...fields], ['symmetry']);
  assert.deepEqual([...bakedTopologyFields(null)], []);
});
