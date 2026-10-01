import { readFile, readdir, writeFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import { isMain } from './exit.mjs';
import { resolve } from 'node:path';
import { BAKED_CONSTANT_IDS, LIVE_TOPOLOGY_FIELD, compileShaderDocument, fixedDerivedBinding } from './shader_workbench.mjs';
import { loadOperatorCatalog } from './pattern_documents.mjs';

const ROOT = fileURLToPath(new URL('../', import.meta.url));
const SECTIONS = ['identity', 'params'];
const GROUPS = { sample: 'source', project: 'projection', warp1: 'outer_warp',
  warp2: 'inner_warp', surface: 'surface', transfer: 'value', cutout: 'value', colorize: 'color' };
const RENAMES = { 'sample.angle-speed': 'angle_rate', 'sample.drift': 'secondary_rate',
  'sample.lattice-shape': 'lattice_shape_blend', 'project.projection-spin-speed': 'spin_rate',
  'project.projection-wander': 'wander', 'colorize.value-opacity-low': 'opacity_low',
  'colorize.value-opacity-high': 'opacity_high' };
const catalog = await loadOperatorCatalog();

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

export function presetAssignments(document, values) {
  const assignments = new Map();
  const topology = topologyIds(document);
  for (const [id, value] of Object.entries(values)) {
    const [label, field] = id.split('.');
    if (field === 'lattice-period') {
      if (!fixedDerivedBinding(document.descriptor, id, values)?.valid)
        throw new Error(`Invalid derived binding: ${id}`);
      continue;
    }
    if (topology.has(id)) continue;
    if (BAKED_CONSTANT_IDS.has(id)) continue;
    let member;
    if (id === 'camera.wander') member = 'projection.camera_wander';
    else if (id === 'sample.edge-width') member = 'value.edge_width';
    else if (/^lens\.mobius-[abcd]-(re|im)$/.test(id))
      member = `lens.mobius.${field.split('-').slice(1).join('.')}`;
    else if (GROUPS[label]) member = `${GROUPS[label]}.${RENAMES[id] ?? field.replaceAll('-', '_')}`;
    else throw new Error(`Unmapped composed parameter: ${id}`);
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
  if (edge.path_policy !== 'parallel') throw new Error(`Unmapped path policy: ${edge.path_policy}`);
  const easing = EASINGS[edge.easing];
  if (!easing) throw new Error(`Unmapped easing: ${edge.easing}`);
  return `Segue::Preset::Lerp{${edge.duration}, ${easing}}`;
}

export function generatedSections(compiled) {
  if (compiled.status !== 'VALID') throw new Error(JSON.stringify(compiled.diagnostics));
  const document = compiled.document;
  const bank = document.preset_bank;
  const order = bank.choreography.generated_order;
  const presets = order.map((id) => bank.presets.find((preset) => preset.preset_id === id));
  const dwell = order.map((id) => bank.choreography.dwell[id]);
  if (!dwell.every((value) => value === dwell[0])) throw new Error('Composed dwell must be uniform');
  const spin = presets.map((preset) => preset.values['camera.spin-speed']);
  if (!spin.every((value) => value === spin[0])) throw new Error('Composed camera spin must be uniform');
  for (const parameter of document.descriptor.parameters) {
    if (parameter.storage !== 'enum8' || parameter.id.endsWith('.palette-mapping')) continue;
    const values = presets.map((preset) => preset.values[parameter.id]);
    if (!values.every((value) => value === values[0]))
      throw new Error(`Composed topology must be uniform: ${parameter.id}`);
  }
  const identity = [
    `  static constexpr std::string_view EFFECT_ID = ${JSON.stringify(document.effect_id)};`,
    `  static constexpr std::string_view DESCRIPTOR_DIGEST = "${compiled.descriptor_digest}";`,
    `  static constexpr std::string_view PRESET_BANK_DIGEST = "${compiled.preset_bank_digest}";`,
    `  static constexpr std::array<std::string_view, ${order.length}> PRESET_IDS{`,
    ...order.map((id, i) => `      ${JSON.stringify(id)}${i + 1 < order.length ? ',' : ''}`),
    '  };',
    `  static constexpr uint16_t PRESET_DWELL_FRAMES = ${dwell[0]};`,
    ...(spin[0] === undefined ? [] : [`  static constexpr float CAMERA_SPIN_RATE = ${floatLiteral(spin[0])};`]),
  ].join('\n');
  const initial = presetAssignments(document, presets[0].values);
  const assignment = ([member, literal]) => `    value.${member} = ${literal};`;
  const params = [
    '  static constexpr Params initial_params() {',
    '    Params value;', ...[...initial].map(assignment), '    return value;', '  }',
  ];
  if (presets.length > 1) {
    const departures = order.map((id) => departureLiteral(bank, id, document.descriptor.path_policies));
    const uniform = departures.every((departure) => departure === departures[0]);
    params.push('', '  /** @brief The preset at index in PRESET_IDS and how it departs. */',
      '  HS_COLD_MEMBER static constexpr PresetEntry<Params> preset(size_t index) {', '    Params value = initial_params();');
    if (!uniform) params.push(`    Segue::Preset::Departure segue = ${departures[0]};`);
    for (let i = 1; i < presets.length; ++i) {
      params.push(`    if (index == ${i}) {`);
      for (const entry of presetAssignments(document, presets[i].values))
        if (initial.get(entry[0]) !== entry[1]) params.push(`  ${assignment(entry)}`);
      if (!uniform && departures[i] !== departures[0]) params.push(`      segue = ${departures[i]};`);
      params.push('    }');
    }
    params.push(`    return {value, ${uniform ? departures[0] : 'segue'}};`, '  }');
  }
  return { identity, params: params.join('\n') };
}

export function updateHeader(source, compiled) {
  const sections = generatedSections(compiled);
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
  let count = 0;
  for (const file of await readdir(resolve(root, 'effects'))) {
    if (!file.endsWith('.h')) continue;
    const path = resolve(root, 'effects', file);
    const source = (await readFile(path, 'utf8')).replaceAll('\r\n', '\n');
    const id = /static constexpr std::string_view EFFECT_ID = "([^"]+)"/.exec(source)?.[1];
    if (!id || !source.includes('DESCRIPTOR_DIGEST')) continue;
    const input = await readFile(resolve(root, 'patterns', `${id.replaceAll('-', '_')}.shader.json`), 'utf8');
    const output = updateHeader(source, compileShaderDocument(input, { catalog }));
    if (check && source !== output) throw new Error(`${file}: regenerate composed presets`);
    if (!check) await writeFile(path, output);
    ++count;
  }
  if (count === 0) throw new Error('No composed presets found');
  return count;
}

if (isMain(import.meta.url)) {
  if (process.argv.slice(2).some((arg) => arg !== '--check')) throw new Error('Usage: generate_composed_presets.mjs [--check]');
  const check = process.argv.includes('--check');
  console.log(`${check ? 'Verified' : 'Regenerated'} ${await generate({ check })} composed preset headers.`);
}
