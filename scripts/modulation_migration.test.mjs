import { test } from 'node:test';
import assert from 'node:assert/strict';
import { declarationFromCatalogField, migrateLegacyShaderDocument } from './shader_workbench.mjs';

const enumField = (id, values, value = values[0]) => ({ id, topology: true, values, default: value });
const amount = { id: 'modulation-amount', min: -4, max: 4, default: 0, curve: 'lerp' };
const modulation = { id: 'color.modulate.v1', params: [amount,
  enumField('modulation-enabled', ['disabled', 'enabled'], 'enabled'),
  enumField('modulation-source', ['noise', 'directional-sine', 'sphere-interference', 'sphere-smooth-grid', 'displacement']),
  enumField('modulation-evaluation', ['direct', 'baked-cube'], 'baked-cube'),
  enumField('modulation-mechanic', ['hue-rotation', 'positive-lightness-lift', 'opacity']),
  { id: 'pattern-speed', min: -0.001, max: 0.001, default: 0, curve: 'lerp' },
] };
const catalog = { operators: [{ id: 'colorize.generated-palette.v4', params: [] }, modulation] };

function legacy() {
  const fields = [
    { ...amount, id: 'hue-shift-amount', default: -0.75 },
    enumField('hue-shift-mode', ['none', 'path-length']),
    ...modulation.params.filter((field) => field.id !== 'modulation-amount'),
    { ...amount, default: 0.4 },
  ];
  const parameters = fields.map((field) => declarationFromCatalogField('ink', field, 'colorize.generated-palette.v3'));
  const values = Object.fromEntries(parameters.map((parameter) => [parameter.id, parameter.default]));
  return { descriptor: {
    chain: [{ label: 'displacement', operator: 'warp.test.v1' },
      { label: 'modulation', operator: 'field.test.v1' },
      { label: 'ink', operator: 'colorize.generated-palette.v3' }],
    parameters,
    path_policies: [{ id: 'serial', kind: 'STAGGERED_ORDERED', groups: parameters.map((parameter) => parameter.id) }],
    serialization: { schema_version: 1, fields: parameters.map((parameter) => parameter.id) },
  }, preset_bank: { presets: [{ preset_id: 'a', values }, { preset_id: 'b', values: {
    ...values, 'ink.hue-shift-mode': 'path-length', 'ink.hue-shift-amount': 1.5,
    'ink.modulation-enabled': 'disabled', 'ink.pattern-speed': -0.0003,
  } }] } };
}

test('legacy colorizers preserve independent path and spatial controls in stack order', () => {
  const input = legacy();
  const before = structuredClone(input);
  const result = migrateLegacyShaderDocument(input, catalog);
  assert.deepEqual(input, before);
  assert.deepEqual(result.descriptor.chain.slice(-3), [
    { label: 'ink', operator: 'colorize.generated-palette.v4' },
    { label: 'displacement-2', operator: 'color.modulate.v1' },
    { label: 'modulation-2', operator: 'color.modulate.v1' },
  ]);
  const [a, b] = result.preset_bank.presets.map((preset) => preset.values);
  assert.equal(a['displacement-2.modulation-enabled'], 'disabled');
  assert.equal(b['displacement-2.modulation-enabled'], 'enabled');
  assert.equal(b['displacement-2.modulation-source'], 'displacement');
  assert.equal(b['displacement-2.modulation-evaluation'], 'direct');
  assert.equal(b['displacement-2.modulation-amount'], 1.5);
  assert.equal(b['modulation-2.modulation-enabled'], 'disabled');
  assert.equal(b['modulation-2.pattern-speed'], -0.0003);
  assert.equal(a['modulation-2.modulation-amount'], 0.4);
  assert.equal(result.descriptor.parameters.find((parameter) => parameter.id === 'displacement-2.modulation-amount').default, -0.75);
  const ids = result.descriptor.parameters.map((parameter) => parameter.id);
  assert.deepEqual(new Set(result.descriptor.path_policies[0].groups), new Set(ids));
  assert.deepEqual(result.descriptor.serialization.fields, ids);
  assert.ok(!JSON.stringify(result).includes('hue-shift'));
  assert.strictEqual(migrateLegacyShaderDocument(result, catalog), result);
});

test('legacy colorizers with implicit controls still retain disabled displacement layers', () => {
  const input = legacy();
  input.descriptor.parameters = [];
  input.descriptor.serialization.fields = [];
  input.descriptor.path_policies = [{ id: 'parallel', kind: 'PARALLEL' }];
  input.preset_bank.presets = [{ preset_id: 'a', values: {} }];
  const result = migrateLegacyShaderDocument(input, catalog);
  assert.equal(result.descriptor.chain.filter((entry) => entry.operator === modulation.id).length, 2);
  assert.equal(result.preset_bank.presets[0].values['displacement-2.modulation-enabled'], 'disabled');
  assert.equal(result.preset_bank.presets[0].values['displacement-2.modulation-amount'], 0);
});

test('standalone legacy modulation stages retain their labels and ordered controls', () => {
  const input = legacy();
  input.descriptor.chain = [
    { label: 'ink', operator: 'colorize.generated-palette.v4' },
    { label: 'hue', operator: 'modulate.hue-shift.v1' },
    { label: 'field', operator: 'modulate.field.v1' },
  ];
  input.descriptor.parameters = [
    declarationFromCatalogField('hue', { ...amount, id: 'hue-shift-amount', default: -0.75 }, 'modulate.hue-shift.v1'),
    declarationFromCatalogField('hue', enumField('hue-shift-mode', ['none', 'path-length'], 'path-length'), 'modulate.hue-shift.v1'),
    ...modulation.params.map((field) => declarationFromCatalogField('field', field, 'modulate.field.v1')),
  ];
  const values = Object.fromEntries(input.descriptor.parameters.map((parameter) => [parameter.id, parameter.default]));
  input.preset_bank.presets = [{ preset_id: 'a', values }, { preset_id: 'b', values: {
    ...values, 'hue.hue-shift-mode': 'none', 'field.modulation-amount': 0.7,
  } }];
  input.descriptor.path_policies[0].groups = Object.keys(values);
  input.descriptor.serialization.fields = Object.keys(values);
  const result = migrateLegacyShaderDocument(input, catalog);
  assert.deepEqual(result.descriptor.chain, [
    { label: 'ink', operator: 'colorize.generated-palette.v4' },
    { label: 'hue', operator: 'color.modulate.v1' },
    { label: 'field', operator: 'color.modulate.v1' },
  ]);
  assert.equal(result.preset_bank.presets[0].values['hue.modulation-amount'], -0.75);
  assert.equal(result.preset_bank.presets[0].values['hue.modulation-enabled'], 'enabled');
  assert.equal(result.preset_bank.presets[1].values['hue.modulation-enabled'], 'disabled');
  assert.equal(result.preset_bank.presets[1].values['field.modulation-amount'], 0.7);
  assert.equal(result.preset_bank.presets[0].values['hue.modulation-source'], 'displacement');
  assert.equal(result.preset_bank.presets[0].values['hue.modulation-evaluation'], 'direct');
  assert.deepEqual(new Set(result.descriptor.path_policies[0].groups),
    new Set(result.descriptor.parameters.map((parameter) => parameter.id)));
  assert.strictEqual(migrateLegacyShaderDocument(result, catalog), result);
});

test('standalone legacy hue modulation keeps its enabled default when mode is implicit', () => {
  const input = legacy();
  input.descriptor.chain = [{ label: 'hue', operator: 'modulate.hue-shift.v1' }];
  input.descriptor.parameters = [];
  input.descriptor.serialization.fields = [];
  input.descriptor.path_policies = [{ id: 'parallel', kind: 'PARALLEL' }];
  input.preset_bank.presets = [{ preset_id: 'a', values: {} }];
  const result = migrateLegacyShaderDocument(input, catalog);
  assert.equal(result.preset_bank.presets[0].values['hue.modulation-enabled'], 'enabled');
  assert.equal(result.preset_bank.presets[0].values['hue.modulation-amount'], 0);
  assert.equal(result.descriptor.chain.length, 1);
});
