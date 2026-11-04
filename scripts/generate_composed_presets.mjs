import { readFile, readdir, writeFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import { exitAfterStderr, isMain } from './exit.mjs';
import { resolve } from 'node:path';
import { LIVE_TOPOLOGY_FIELD, compileShaderDocument, fixedDerivedBinding } from './shader_workbench.mjs';
import { loadOperatorCatalog } from './pattern_documents.mjs';
import { loadEffectIdentities } from './effect_roster.mjs';

const ROOT = fileURLToPath(new URL('../', import.meta.url));
const SECTIONS = ['identity', 'params'];
const GROUPS = { sample: 'source', project: 'projection', warp1: 'outer_warp',
  warp2: 'inner_warp', warp: 'outer_warp', displace: 'displace', surface: 'surface', transfer: 'value', cutout: 'value', colorize: 'color' };
const MODULATION_MEMBERS = {
  'modulation-amount': 'mechanic.amount',
  'pattern-scale': 'noise.scale',
  'pattern-speed': 'noise.speed',
  'pattern-contrast': 'noise.contrast',
};
const MODULATION_ENUMS = {
  'modulation-source': ['source', 'Pattern::SourceKind'],
  'modulation-evaluation': ['evaluation', 'Pattern::Evaluation'],
  'modulation-mechanic': ['mechanic.kind', 'ColorModulation::Mechanic'],
  'pattern-noise-basis': ['noise.basis', 'Pattern::NoiseBasis'],
  'pattern-noise-variant': ['noise.variant', 'Pattern::NoiseVariant'],
};
const catalog = await loadOperatorCatalog();
const identities = await loadEffectIdentities();

const LIVE_SURFACE_OPERATORS = new Set(['sphere.displace.symmetric.v2', 'sphere.displace.rossby-haurwitz.v1', 'sphere.displace.ripple.v2']);

const liveTopology = (document, id) => {
  const [label, field] = id.split('.');
  const slot = document.descriptor.chain.find(slot => slot.label === label);
  return LIVE_SURFACE_OPERATORS.has(slot?.operator) &&
    catalog.operators.find(operator => operator.id === slot.operator).params
      .find(parameter => parameter.id === field)?.topology;
};

const authoredLiveTopology = (slot, field) => {
  if (slot.operator.startsWith('sample.')) return ['basis', 'kind', 'weight-mode', 'coverage-mode'].includes(field);
  return slot.operator !== 'sphere.rotate.v2';
};

const topologyIds = (document) => new Set(document.descriptor.chain.flatMap((slot) =>
  catalog.operators.find((operator) => operator.id === slot.operator).params
    .filter((field) => field.topology && field.id !== LIVE_TOPOLOGY_FIELD)
    .map((field) => `${slot.label}.${field.id}`)));

export function floatLiteral(value) {
  const rounded = Math.fround(value);
  if (!Number.isFinite(rounded)) throw new Error('Nonfinite preset value');
  if (Object.is(rounded, -0)) return '-0.0f';
  let text;
  for (let digits = 1; digits <= 9; ++digits) {
    text = Number(rounded.toPrecision(digits)).toString();
    if (Math.fround(Number(text)) === rounded) break;
  }
  return `${/[.e]/.test(text) ? text : `${text}.0`}f`;
}

const groupFor = (slot, authoredTopology = false) => {
  if (slot.resource_key) return slot.resource_key;
  if (slot.operator === 'sphere.rotate.v2' && authoredTopology) return slot.label;
  if (slot.operator.startsWith('sample.')) return 'source';
  if (slot.operator.startsWith('project.') || slot.operator === 'sphere.rotate.v2') return 'projection';
  if (slot.operator.startsWith('colorize.')) return 'color';
  return GROUPS[slot.label] ?? slot.label;
};

const scalarMember = (slot, field, authoredTopology = false) => {
  if (slot.operator === 'sphere.rotate.v2')
    return authoredTopology ? field.replaceAll('-', '_') : field === 'wander' ? 'camera_wander' : null;
  if (slot.operator.startsWith('sample.')) {
    const sourceNames = { 'angle-speed': 'angle_rate', drift: 'secondary_rate',
      'lattice-shape': 'lattice_shape_blend', 'noise-speed': 'noise_time_rate' };
    return sourceNames[field] ?? field.replaceAll('-', '_');
  }
  if (slot.operator.startsWith('project.'))
    return ({ 'projection-spin-speed': 'spin_rate', 'projection-wander': 'wander' })[field] ?? field.replaceAll('-', '_');
  if (slot.operator === 'sphere.lens.twist.v2' && field === 'twist-distribution') return 'distribution';
  if (slot.operator === 'sphere.displace.potential.v1' && /^term-[1-4]-/.test(field)) {
    const [, index, ...parts] = field.split('-');
    return `terms[${Number(index) - 1}].${parts.join('_')}`;
  }
  return ({ 'value-opacity-low': 'opacity_low', 'value-opacity-high': 'opacity_high' })[field] ?? field.replaceAll('-', '_');
};

/** @brief Catalog choices pinned by a concrete pipeline's policy types. */
export function topologySelections(document, values) {
  const selections = new Map();
  for (const slot of document.descriptor.chain) {
    const operator = catalog.operators.find(operator => operator.id === slot.operator);
    if (!operator) throw new Error(`Unknown composed operator: ${slot.operator}`);
    for (const field of operator.params.filter(field => field.topology)) {
      const id = `${slot.label}.${field.id}`;
      if (!(id in values)) continue;
      if (!field.values.includes(values[id])) throw new Error(`Invalid topology selector: ${id}`);
      selections.set(id, values[id]);
    }
  }
  return selections;
}

export function presetAssignments(document, values, { authoredTopology = false } = {}) {
  const assignments = new Map();
  const topology = topologyIds(document);
  const slots = new Map(document.descriptor.chain.map(slot => [slot.label, slot]));
  for (const [id, value] of Object.entries(values)) {
    const [label, field] = id.split('.');
    const slot = slots.get(label);
    const definition = catalog.operators.find(operator => operator.id === slot?.operator)?.params
      .find(parameter => parameter.id === field);
    if (!slot || !definition) throw new Error(`Unmapped composed parameter: ${id}`);
    const group = groupFor(slot, authoredTopology);
    if (field === 'lattice-period') {
      if (!['warp.affine.v3', 'warp.affine-periodic.v1'].includes(slot.operator)) throw new Error(`Invalid derived binding: ${id}`);
      if (!authoredTopology && slot.operator === 'warp.affine.v3' && fixedDerivedBinding(document.descriptor, id, values)?.valid) continue;
      if (typeof value !== 'number' || !Number.isFinite(value) || value <= 0)
        throw new Error(`Invalid affine lattice period: ${id}`);
      assignments.set(`${group}.lattice_period`, floatLiteral(value));
      continue;
    }
    if (slot.operator === 'color.modulate.v1') {
      if (field === 'modulation-enabled') {
        if (!['disabled', 'enabled'].includes(value)) throw new Error(`Invalid modulation state: ${value}`);
        assignments.set(`${group}.modulation.enabled`, String(value === 'enabled'));
        continue;
      }
      if (MODULATION_ENUMS[field]) {
        const [path, type] = MODULATION_ENUMS[field];
        if (!definition.values.includes(value)) throw new Error(`Invalid modulation selector: ${id}`);
        assignments.set(`${group}.modulation.${path}`, `${type}::${value.replaceAll('-', '_').toUpperCase()}`);
        continue;
      }
      if (MODULATION_MEMBERS[field]) {
        assignments.set(`${group}.modulation.${MODULATION_MEMBERS[field]}`, floatLiteral(value));
        continue;
      }
      if (/^(sine|interference|grid)-/.test(field)) {
        const [source, ...parts] = field.split('-');
        const path = parts.join('_').replace(/^axis_b_([xyz])$/, 'axis_b.$1').replace(/^axis_([xyz])$/, 'axis.$1');
        assignments.set(`${group}.modulation.${source}.${path}`, floatLiteral(value));
        continue;
      }
    }
    if (topology.has(id)) {
      if (!definition.values.includes(value)) throw new Error(`Invalid topology selector: ${id}`);
      if (liveTopology(document, id) || (authoredTopology && authoredLiveTopology(slot, field))) {
        const member = field === 'brightness-envelope' ? 'brightness_envelope' : scalarMember(slot, field);
        const literal = field === 'brightness-envelope'
          ? `Pullback::Color::BrightnessEnvelope::${value.replaceAll('-', '_').toUpperCase()}`
          : String(definition.values.indexOf(value));
        assignments.set(`${group}.${member}`, literal);
      }
      continue;
    }
    if (!authoredTopology && slot.operator === 'sphere.rotate.v2' && field === 'spin-speed') continue;
    if (slot.operator.startsWith('sample.') && field === 'edge-width') {
      if (authoredTopology) {
        assignments.set(`${group}.edge_width`, floatLiteral(value));
        continue;
      }
      if (!fixedDerivedBinding(document.descriptor, id, values))
        assignments.set('value.edge_width', floatLiteral(value));
      continue;
    }
    const member = `${group}.${scalarMember(slot, field, authoredTopology)}`;
    if (field === 'palette-mapping') {
      if (!['cup', 'bell', 'linear', 'reverse'].includes(value)) throw new Error(`Invalid mapping: ${value}`);
      assignments.set(member, `Pullback::Color::PaletteMapping::${value.toUpperCase()}`);
    } else {
      if (typeof value !== 'number') throw new Error(`Unmapped composed value: ${id}`);
      assignments.set(member, floatLiteral(value));
    }
  }
  return assignments;
}

const EASINGS = { LINEAR: 'math::ease_linear', EASE_IN_OUT_SIN: 'math::ease_in_out_sin' };

/** The C++ departure of a preset: the document edge leaving it, or a snap. */
export function departureLiteral(bank, presetId, policies) {
  const edges = bank.edges.filter((edge) => edge.from === presetId);
  if (edges.length === 0) {
    if (bank.absent_edge_fallback.automatic !== 'SNAP')
      throw new Error('Absent automatic edge must fall back to SNAP');
    return 'Segue::Preset::Snap{}';
  }
  if (edges.length > 1) throw new Error(`Preset ${presetId} departs by several edges`);
  const [edge] = edges;
  const order = bank.choreography.generated_order;
  const index = order.indexOf(presetId);
  if (index < 0 || edge.to !== order[(index + 1) % order.length])
    throw new Error(`Preset ${presetId} edge does not target its successor`);
  if (policies.find(policy => policy.id === edge.path_policy)?.kind !== 'PARALLEL')
    throw new Error(`Unmapped path policy kind: ${edge.path_policy}`);
  const easing = EASINGS[edge.easing];
  if (!easing) throw new Error(`Unmapped easing: ${edge.easing}`);
  return `Segue::Preset::Lerp{${edge.duration}, ${easing}}`;
}

export function generatedSections(compiled, metadata, { authoredTopology = false } = {}) {
  if (compiled.status !== 'VALID') throw new Error(JSON.stringify(compiled.diagnostics));
  const document = compiled.document;
  const bank = document.preset_bank;
  const order = bank.choreography.generated_order;
  const presets = order.map((id) => bank.presets.find((preset) => preset.preset_id === id));
  if (!bank.choreography.dwell)
    throw new Error('Composed documents must declare choreography.dwell');
  const dwell = order.map((id) => bank.choreography.dwell[id]);
  if (!dwell.every((value) => value === dwell[0])) throw new Error('Composed dwell must be uniform');
  const camera = document.descriptor.chain.find(slot => slot.operator === 'sphere.rotate.v2');
  const spin = presets.map(preset => camera ? preset.values[`${camera.label}.spin-speed`] : undefined);
  if (!authoredTopology && !spin.every((value) => value === spin[0])) throw new Error('Composed camera spin must be uniform');
  for (const parameter of document.descriptor.parameters) {
    const [label, field] = parameter.id.split('.');
    const slot = document.descriptor.chain.find(slot => slot.label === label);
    if (parameter.storage !== 'enum8' || field === LIVE_TOPOLOGY_FIELD ||
      liveTopology(document, parameter.id) || slot.operator === 'color.modulate.v1' ||
      (authoredTopology && authoredLiveTopology(slot, field))) continue;
    const values = presets.map((preset) => preset.values[parameter.id]);
    if (!values.every((value) => value === values[0]))
      throw new Error(`Composed topology must be uniform: ${parameter.id}`);
  }
  metadata ??= identities.find(row => row.id === document.effect_id);
  if (!metadata) throw new Error(`Unregistered composed effect ID: ${document.effect_id}`);
  const resources = new Map();
  if (authoredTopology)
    for (const slot of document.descriptor.chain) {
      const resource = groupFor(slot, true);
      if (resources.has(resource) && resources.get(resource) !== slot.label)
        throw new Error(`Several authored instances share resource "${resource}"`);
      resources.set(resource, slot.label);
    }
  const identity = [
    '  /// Machine/display identity and compatibility keys.',
    `  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::${metadata.name}::VALUE;`,
    "  /// Canonical machine ID; equals the shader document's effect_id.",
    '  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;',
    "  /// SHA-256 hex of the shader document's parameter descriptor.",
    `  static constexpr std::string_view DESCRIPTOR_DIGEST = "${compiled.descriptor_digest}";`,
    "  /// SHA-256 hex of the shader document's preset bank.",
    `  static constexpr std::string_view PRESET_BANK_DIGEST = "${compiled.preset_bank_digest}";`,
    '  /// Preset identities, indexed by preset number.',
    `  static constexpr std::array<std::string_view, ${order.length}> PRESET_IDS{`,
    ...order.map((id, i) => `      ${JSON.stringify(id)}${i + 1 < order.length ? ',' : ''}`),
    '  };',
    '  /// Frames each preset holds before advancing.',
    `  static constexpr uint16_t PRESET_DWELL_FRAMES = ${dwell[0]};`,
    ...(!authoredTopology ? [] : [
      '  /** @brief Authored instance label of a native resource.',
      '   * @param resource Native resource key.',
      '   * @return Descriptor instance label. */',
      '  static constexpr std::string_view authored_instance(std::string_view resource) {',
      ...[...resources].filter(([resource, label]) => resource !== label)
        .map(([resource, label]) => `    if (resource == ${JSON.stringify(resource)}) return ${JSON.stringify(label)};`),
      '    return resource;',
      '  }',
    ]),
    ...(authoredTopology || spin[0] === undefined ? [] : [
      '  /// Outer camera spin about +Y, in radians per frame.',
      `  static constexpr float CAMERA_SPIN_RATE = ${floatLiteral(spin[0])};`,
    ]),
  ].join('\n');
  const initial = presetAssignments(document, presets[0].values, { authoredTopology });
  const assignment = ([member, literal]) => {
    const [key, ...path] = member.split('.');
    return `    value.template get<"${key}">().${path.join('.')} = ${literal};`;
  };
  const params = [
    '  /**',
    '   * @brief Parameters of the first preset; every preset varies from them.',
    '   * @return The preset-0 `Params`.',
    '   */',
    '  static constexpr Params initial_params() {',
    '    Params value;', ...[...initial].map(assignment), '    return value;', '  }',
  ];
  if (presets.length > 1) {
    const departures = order.map((id) => departureLiteral(bank, id, document.descriptor.path_policies));
    const uniform = departures.every((departure) => departure === departures[0]);
    params.push('', '  /**',
      '   * @brief The preset at index in PRESET_IDS and how it departs.',
      '   * @param index Preset number; an index past the last yields preset 0.',
      "   * @return The preset's parameters and departure segue.",
      '   */',
      '  HS_COLD_MEMBER static constexpr PresetEntry<Params> preset(size_t index) {', '    Params value = initial_params();');
    if (!uniform) params.push(`    Segue::Preset::Departure segue = ${departures[0]};`);
    for (let i = 1; i < presets.length; ++i) {
      params.push(`    if (index == ${i}) {`);
      for (const entry of presetAssignments(document, presets[i].values, { authoredTopology }))
        if (initial.get(entry[0]) !== entry[1]) params.push(`  ${assignment(entry)}`);
      if (!uniform && departures[i] !== departures[0]) params.push(`      segue = ${departures[i]};`);
      params.push('    }');
    }
    params.push(`    return {value, ${uniform ? departures[0] : 'segue'}};`, '  }');
  }
  return { identity, params: params.join('\n') };
}

export function updateHeader(source, compiled, metadata) {
  const authoredTopology = /AUTHORED_TOPOLOGY\s*=\s*true/u.test(source);
  const sections = generatedSections(compiled, metadata, { authoredTopology });
  for (const section of SECTIONS) {
    const begin = `  // Generated ${section}: scripts/generate_composed_presets.mjs\n  // clang-format off\n`;
    const end = `  // clang-format on\n  // End generated ${section}.`;
    const start = source.indexOf(begin);
    const finish = source.indexOf(end, start);
    if (start < 0 || finish < 0) throw new Error(`Missing generated ${section} markers`);
    source = source.slice(0, start + begin.length) + sections[section] + '\n' + source.slice(finish);
  }
  return source;
}

export async function generate({ check = false, root = ROOT } = {}) {
  const metadata = await loadEffectIdentities(root);
  const composed = [];
  const ids = new Set();
  const stale = [];
  const emit = async (file, output) => {
    const path = resolve(root, file);
    if (!check) return writeFile(path, output);
    let source;
    try {
      source = (await readFile(path, 'utf8')).replaceAll('\r\n', '\n');
    } catch (error) {
      if (error.code !== 'ENOENT') throw error;
    }
    if (source !== output) stale.push(file);
  };
  for (const file of (await readdir(resolve(root, 'effects'))).sort()) {
    if (!file.endsWith('.h')) continue;
    const path = resolve(root, 'effects', file);
    const source = (await readFile(path, 'utf8')).replaceAll('\r\n', '\n');
    if (!source.includes('DESCRIPTOR_DIGEST')) continue;
    const name = /EFFECT_IDENTITY\s*=\s*hs::effect_identities::(\w+)/u.exec(source)?.[1];
    const declaredId = /static constexpr std::string_view EFFECT_ID = "([^"]+)"/u.exec(source)?.[1];
    const record = metadata.find(row => name ? row.name === name : row.id === declaredId);
    if (!record) throw new Error(`Unregistered composed identity: effects/${file}`);
    const id = record.id;
    if (ids.has(id)) throw new Error(`Duplicate composed effect ID: ${id}`);
    ids.add(id);
    const input = await readFile(resolve(root, 'patterns', `${id.replaceAll('-', '_')}.shader.json`), 'utf8');
    const compiled = compileShaderDocument(input, { catalog });
    if (compiled.status !== 'VALID')
      throw new Error(`patterns/${id.replaceAll('-', '_')}.shader.json: ${JSON.stringify(compiled.diagnostics)}`);
    if (compiled.document.effect_id !== id)
      throw new Error(`patterns/${id.replaceAll('-', '_')}.shader.json names effect "${compiled.document.effect_id}", expected "${id}"`);
    const output = updateHeader(source, compiled, record);
    await emit(`effects/${file}`, output);
    composed.push(file);
  }
  if (composed.length === 0) throw new Error('No composed presets found');
  await emit('targets/composed_effect_includes.h', [
    '/*',
    ' * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.',
    ' * Licensed under the PolyForm Noncommercial License 1.0.0',
    ' */',
    '#pragma once', '',
    '/**',
    ' * @file composed_effect_includes.h',
    ' * @brief Includes every composed effect header.',
    ' */', '',
    '// Generated by scripts/generate_composed_presets.mjs.',
    ...composed.map(file => `#include "effects/${file}"`), '',
  ].join('\n'));
  await emit('scripts/composed_effect_roster.mjs', [
    '// Generated by scripts/generate_composed_presets.mjs.',
    'export const SHADER_DOCUMENT_EFFECTS = Object.freeze([',
    ...[...ids].sort().map(id => `  '${id}',`),
    ']);', '',
  ].join('\n'));
  if (stale.length) throw new Error(`Stale composed artifacts:\n${stale.join('\n')}\nRegenerate with: node scripts/generate_composed_presets.mjs`);
  return composed.length;
}

if (isMain(import.meta.url)) {
  if (process.argv.slice(2).some((arg) => arg !== '--check')) {
    console.error('Usage: generate_composed_presets.mjs [--check]');
    await exitAfterStderr(2);
  }
  const check = process.argv.includes('--check');
  console.log(`${check ? 'Verified' : 'Regenerated'} ${await generate({ check })} composed preset headers and rosters.`);
}
