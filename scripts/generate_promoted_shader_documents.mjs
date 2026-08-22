import { readdir, readFile, writeFile } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  compileShaderDocument,
  exportShaderDocumentJson,
} from './shader_workbench.mjs';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');

// --check compiles the specs and compares them against the committed
// documents, exiting non-zero on drift; the default rewrites patterns/.
const argv = process.argv.slice(2);
const CHECK = argv.includes('--check');
const unknown = argv.filter((arg) => arg !== '--check');
if (unknown.length) {
  console.error(`unknown argument: ${unknown[0]}`);
  console.error('usage: generate_promoted_shader_documents.mjs [--check]');
  process.exit(2);
}
const TAU = Math.fround(Math.PI * 2);

// Choreography mirrors the shared ComposedEffect cadence.
const PRESET_DWELL_FRAMES = 600;
const INHERITED_SEGUE_FRAMES = 480;

const defaults = Object.freeze({
  grid: {
    'pattern-freq': 1, speed: 0, complexity: 0, 'pattern-mix': 0,
    drift: 0, 'source-angle-speed': 0,
  },
  twin: {
    'pattern-freq': 1, speed: 0, drift: 0, 'source-angle-speed': 0,
  },
  lattice: {
    'lattice-cell-scale': 1, 'lattice-shape': 0,
    'lattice-softness': 0.05, 'lattice-radius': 0.25,
  },
  projection: {
    'pole-fade': 1, 'projection-spin-speed': 0, 'projection-wander': 0,
    'camera-wander': 0, 'central-meridian': 0,
  },
  mirror: {
    speed: 0, rotation: 0, 'cell-x': 1, 'cell-y': 1,
    'offset-x': 0, 'offset-y': 0,
  },
  wave: { speed: 0, strength: 0, frequency: 1, 'field-angle': 0 },
  vector: { speed: 0, strength: 0, scale: 1, 'vector-angle': 0 },
  affine: {
    speed: 0, 'rotation-rate': 0, 'translation-x': 0, 'translation-y': 0,
    'scale-x': 1, 'scale-y': 1, shear: 0,
  },
  polar: { speed: 0, 'radial-scale': 1, 'radial-phase': 0, 'angular-phase': 0 },
  edge: { 'edge-width': 0.1 },
  iso: { 'iso-level': 0.5, 'iso-width': 0.05 },
  color: {
    'palette-chroma': 0.62, 'palette-mapping': 'linear',
    'mapping-frequency': 1, 'mapping-phase': 0,
    'phase-oscillation-depth': 0, 'phase-oscillation-speed': 0,
    'brightness-depth': 1, 'value-opacity-low': 1, 'value-opacity-high': 1,
    'hue-shift-amount': 0, 'hue-noise-scale': 1, 'hue-noise-speed': 0,
  },
  mobius: {
    'mobius-a-re': Math.SQRT1_2, 'mobius-a-im': 0,
    'mobius-b-re': 0, 'mobius-b-im': 0,
    'mobius-c-re': 0, 'mobius-c-im': 0,
    'mobius-d-re': Math.SQRT1_2, 'mobius-d-im': 0,
  },
});

const prefix = (name, values) => Object.fromEntries(
  Object.entries(values).map(([key, value]) => [`${name}-${key}`, value]),
);

const parameterSpec = (id, value, hue, source) => {
  if (id === 'palette-mapping') return {
    id, binding: 'color.palette-mapping', classification: 'preset',
    storage: 'enum8', unit: 'mapping',
    domain: { values: ['cup', 'bell', 'linear', 'reverse'] },
    interpolation: { kind: 'MIXED_ENUM' }, default: value,
  };
  const angle = id.endsWith('rotation') || id.endsWith('field-angle') ||
    id.endsWith('vector-angle') || id.endsWith('radial-phase') ||
    id.endsWith('angular-phase') || id === 'central-meridian';
  const positive = id.endsWith('cell-x') || id.endsWith('cell-y') ||
    id.endsWith('scale-x') || id.endsWith('scale-y') || id.endsWith('radial-scale') ||
    id === 'lattice-cell-scale' || id === 'lattice-softness' ||
    id === 'lattice-radius' || id.endsWith('frequency') || id === 'hue-noise-scale';
  let domain = { minimum: -30, maximum: 30 };
  let unit = 'ratio';
  if (id === 'pattern-freq')
    domain = { minimum: 0.1, maximum: source === 'grid' ? 64 : 20 };
  else if (id === 'speed') domain = { minimum: 0, maximum: 0.5 };
  else if (id === 'source-angle-speed') domain = { minimum: 0, maximum: Math.fround(0.05) };
  else if (id.endsWith('-speed')) domain = id === 'hue-noise-speed'
    ? { minimum: -0.001, maximum: 0.001 }
    : id === 'projection-spin-speed'
      ? { minimum: 0, maximum: 0.05 }
      : { minimum: -0.02, maximum: 0.02 };
  else if (id === 'drift') domain = { minimum: 0, maximum: 1.25 };
  else if (id === 'complexity') domain = { minimum: 0, maximum: 3 };
  else if (id === 'pattern-mix' || id.endsWith('wander') || id.endsWith('depth') ||
           id.endsWith('opacity-low') || id.endsWith('opacity-high') ||
           id === 'lattice-shape' || id === 'iso-level' || id === 'edge-width')
    domain = { minimum: 0, maximum: 1 };
  else if (id === 'pole-fade') domain = { minimum: 1, maximum: 20 };
  else if (id.endsWith('cell-x') || id.endsWith('cell-y'))
    domain = { minimum: 1 / 64, maximum: 8 };
  else if (id.endsWith('offset-x') || id.endsWith('offset-y'))
    domain = { minimum: -8, maximum: 8 };
  else if (id.endsWith('scale-x') || id.endsWith('scale-y') ||
           id.endsWith('radial-scale')) domain = { minimum: 1 / 64, maximum: 64 };
  else if (id === 'lattice-cell-scale') domain = { minimum: 1 / 64, maximum: 8 };
  else if (id === 'lattice-softness' || id === 'iso-width')
    domain = { minimum: 1 / 1024, maximum: 1 };
  else if (id === 'lattice-radius') domain = { minimum: 1 / 64, maximum: 0.49 };
  else if (id.endsWith('frequency')) domain = { minimum: 0.01, maximum: 32 };
  else if (id.endsWith('rotation-rate')) domain = { minimum: -TAU, maximum: TAU };
  else if (id.endsWith('translation-x') || id.endsWith('translation-y') ||
           id.endsWith('shear')) domain = { minimum: -4, maximum: 4 };
  else if (angle) domain = id === 'central-meridian'
    ? { minimum: 0, maximum: TAU } : { minimum: -TAU, maximum: TAU };
  else if (id === 'palette-chroma') domain = { minimum: 0, maximum: 1 };
  else if (id === 'mapping-frequency') domain = { minimum: 1, maximum: 32 };
  else if (id === 'mapping-phase') domain = { minimum: -1, maximum: 1 };
  else if (id === 'phase-oscillation-speed') domain = { minimum: -0.01, maximum: 0.01 };
  else if (id === 'hue-shift-amount') domain = hue === 'path-length'
    ? { minimum: -4, maximum: 4 } : { minimum: -1, maximum: 1 };
  else if (id === 'hue-noise-scale') domain = { minimum: 1 / 64, maximum: 8 };
  else if (id.startsWith('mobius-')) domain = { minimum: -4, maximum: 4 };
  if (id.includes('speed')) unit = 'turn-per-frame';
  else if (angle || id.endsWith('rotation-rate')) unit = 'radian';
  return {
    id, binding: id.replaceAll('-', '.'), classification: 'preset',
    storage: 'binary32', unit, domain,
    interpolation: angle
      ? { kind: 'SHORTEST_PERIODIC', period: TAU }
      : { kind: positive ? 'LOG_POSITIVE' : 'LINEAR' },
    default: value,
  };
};

const stageGraph = (spec) => {
  const resources = [`${spec.palette}-palette`];
  if (spec.hue === 'noise') resources.push('hue-noise');
  const nodes = [
    { label: 'camera', role: 'outer_camera', operator: 'pullback.outer_camera.v1',
      policy: { base_orientation: 'negative-z-to-negative-y' } },
    { label: 'surface', role: 'surface_project', operator: 'pullback.surface_project.v1',
      policy: { pre_lens_surface: 'identity', lens: spec.lens,
        post_lens_surface: 'identity', projection: spec.projection } },
    { label: 'warp', role: 'planar_warp', operator: 'pullback.planar_warp.v1',
      policy: { sequence: [spec.outer, spec.inner] } },
    { label: 'source', role: 'source', operator: 'pullback.source.v1',
      policy: { source: spec.source } },
    { label: 'material', role: 'material', operator: 'pullback.material.v1',
      policy: { weight: 'projection', transfer: spec.transfer,
        coverage: spec.coverage } },
    { label: 'color', role: 'color', operator: 'pullback.color.v1',
      policy: { color: 'generated-palette', hue_mode: spec.hue,
        brightness_envelope: spec.brightness }, resources },
  ];
  const edges = nodes.slice(0, -1).map((node, index) => ({
    from: node.label, to: nodes[index + 1].label,
  }));
  return { nodes, edges };
};

const sourceValues = (kind) => ({ ...defaults[kind] });
const warpValues = (kind, slot) => kind === 'none' ? {} : prefix(slot, defaults[kind]);

const baseValues = (spec) => {
  const values = {
    ...sourceValues(spec.sourceKey),
    'pole-fade': defaults.projection['pole-fade'],
    'camera-wander': defaults.projection['camera-wander'],
    ...warpValues(spec.outerKey, 'outer'),
    ...warpValues(spec.innerKey, 'inner'),
    ...(spec.valueKey ? defaults[spec.valueKey] : {}),
    'palette-chroma': defaults.color['palette-chroma'],
    'palette-mapping': defaults.color['palette-mapping'],
    'mapping-frequency': defaults.color['mapping-frequency'],
    'mapping-phase': defaults.color['mapping-phase'],
    'phase-oscillation-depth': defaults.color['phase-oscillation-depth'],
    'phase-oscillation-speed': defaults.color['phase-oscillation-speed'],
    'value-opacity-low': defaults.color['value-opacity-low'],
    'value-opacity-high': defaults.color['value-opacity-high'],
    'hue-shift-amount': defaults.color['hue-shift-amount'],
  };
  if (spec.animatedProjection) {
    values['projection-spin-speed'] = 0;
    values['projection-wander'] = 0;
  }
  if (spec.centralMeridian) values['central-meridian'] = 0;
  if (spec.hue === 'noise') {
    values['hue-noise-scale'] = 1;
    values['hue-noise-speed'] = 0;
  }
  if (spec.brightness !== 'none') values['brightness-depth'] = 1;
  if (spec.lens === 'mobius') Object.assign(values, defaults.mobius);
  return values;
};

const bank = (spec, base) => {
  const presets = spec.presets.map((preset) => ({
    preset_id: preset.id,
    display_name: preset.name,
    values: { ...base, ...preset.values },
  }));
  const edges = presets.length < 2 ? []
    : presets.map((preset, index) => ({
      from: preset.preset_id,
      to: presets[(index + 1) % presets.length].preset_id,
      path_policy: 'parallel', easing: 'EASE_IN_OUT_SIN',
      duration: INHERITED_SEGUE_FRAMES,
    }));
  return {
    schema_version: 1, presets, edges,
    absent_edge_fallback: {
      manual: 'SNAP', automatic: 'REJECT', synchronized: 'SNAP',
      restore: 'SNAP', authoring: 'SNAP',
    },
    choreography: {
      generated_order: presets.map((preset) => preset.preset_id),
      dwell: Object.fromEntries(presets.map(
        (preset) => [preset.preset_id, PRESET_DWELL_FRAMES])),
    },
  };
};

const documentFor = (spec) => {
  const values = baseValues(spec);
  const parameters = Object.entries(values)
    .map(([id, value]) => parameterSpec(id, value, spec.hue, spec.source));
  const resources = [{ id: `${spec.palette}-palette`, kind: `generated-${spec.palette}-palette`,
    settings: { hue_step: 159 } }];
  if (spec.hue === 'noise') resources.push({
    id: 'hue-noise', kind: 'simplex-sphere-scalar', settings: { seed: 6047 },
  });
  if (spec.outerKey === 'vector') resources.push({
    id: 'outer-warp-noise', kind: 'simplex-planar-vector', settings: { seed: 1337 },
  });
  return {
    schema_version: 1,
    catalog_version: 1,
    document_id: `${spec.id}-v1`,
    effect_id: spec.id,
    descriptor: {
      graph: stageGraph(spec),
      parameters,
      path_policies: [{ id: 'parallel', kind: 'PARALLEL' }],
      clocks: [
        { id: 'palette-cycle', kind: 'generated-palette-cycle', settings: { fade: 600 } },
        { id: 'transition-progress', kind: 'completed-frame', settings: { pause: 'parameter-animation' } },
      ],
      preparation: [
        { id: 'spatial-frames', kind: 'prepare-orientation-conjugates' },
        { id: 'hue-luts', kind: 'prepare-hue-luts' },
        { id: 'warp-carriers', kind: 'prepare-fixed-warps' },
      ],
      resources,
      serialization: { schema_version: 1, fields: parameters.map((parameter) => parameter.id) },
      approximation: [{ id: 'hue-luts', kind: 'registered-oracle', settings: { oracle: 'hue-luts-v1' } }],
      handoff: { policy: 'reset' },
    },
    preset_bank: bank(spec, values),
    effect_metadata: { display_name: spec.display, description: spec.description },
  };
};

// The twelve promoted documents this generator owns. patterns/ holds three
// more — curl_lattice, facet_grid and prism_spiral — hand-authored from
// workbench snapshots, which no rerun writes.
const effects = [
  {
    id: 'signal-weave', display: 'Signal Weave', source: 'grid', sourceKey: 'grid',
    projection: 'stereographic', lens: 'glitch', outer: 'wave-shear', outerKey: 'wave',
    inner: 'identity', innerKey: 'none', transfer: 'linear', coverage: 'projection-squared',
    palette: 'triadic', hue: 'noise', brightness: 'none', animatedProjection: true,
    description: 'Glitch-folded grids pulled through an animated wave shear.',
    presets: [
      { id: 'signal-weave', name: 'Signal Weave', values: {
        'pattern-freq': 4.439, speed: 0.245, complexity: 0.5,
        'camera-wander': 0.8, 'outer-strength': 0.5, 'outer-speed': 0.015625,
        'hue-shift-amount': 0.292, 'hue-noise-scale': 0.6304219,
        'palette-chroma': 0.788, 'palette-mapping': 'cup',
      } },
      { id: 'signal-weave-2', name: 'Signal Weave 2', values: {
        'pattern-freq': 3.1447, speed: 0.245, complexity: 0.5,
        'camera-wander': 0.8, 'outer-strength': 2.72, 'outer-speed': 0.00690625,
        'hue-shift-amount': 0.292, 'hue-noise-scale': 0.6304219,
        'palette-chroma': 0.788, 'palette-mapping': 'cup',
      } },
      { id: 'signal-weave-3', name: 'Signal Weave 3', values: {
        'pattern-freq': 7.5227, speed: 0.245, complexity: 1.698,
        'camera-wander': 0.8, 'outer-strength': 0, 'outer-speed': 0.00690625,
        'hue-shift-amount': 0.292, 'hue-noise-scale': 0.6304219,
        'palette-chroma': 0.788, 'palette-mapping': 'cup',
      } },
      { id: 'signal-weave-4', name: 'Signal Weave 4', values: {
        'pattern-freq': 8.8162, speed: 0.245, complexity: 1.698,
        'camera-wander': 0.8, 'outer-strength': 1.376, 'outer-speed': 0.00559375,
        'hue-shift-amount': 0.292, 'hue-noise-scale': 0.6304219,
        'palette-chroma': 0.788, 'palette-mapping': 'cup',
      } },
    ],
  },
  {
    id: 'kaleido-wave', display: 'Kaleido Wave', source: 'twin-wave', sourceKey: 'twin',
    projection: 'stereographic', lens: 'kaleidoscope', outer: 'identity', outerKey: 'none',
    inner: 'mirror-tile', innerKey: 'mirror', transfer: 'linear', coverage: 'projection-squared',
    palette: 'triadic', hue: 'noise', brightness: 'none', animatedProjection: true,
    description: 'A drifting twin wave reflected through a kaleidoscope.',
    presets: [{ id: 'twin-wave', name: 'Twin Wave', values: {
      'pattern-freq': 4.9755, speed: 0.125, drift: 0.8, 'source-angle-speed': 0.05,
      'pole-fade': 4.971, 'projection-wander': 1, 'camera-wander': 1,
      'hue-shift-amount': 0.27, 'hue-noise-scale': 2.2033439,
      'hue-noise-speed': -0.00040800002, 'palette-chroma': 0.361,
    } }],
  },
  {
    id: 'alien-ocean', display: 'Alien Ocean', source: 'grid', sourceKey: 'grid',
    projection: 'gnomonic-folded', lens: 'kaleidoscope', outer: 'mirror-tile', outerKey: 'mirror',
    inner: 'identity', innerKey: 'none', transfer: 'linear', coverage: 'edge-fade', valueKey: 'edge',
    palette: 'triadic', hue: 'noise', brightness: 'none', animatedProjection: false,
    description: 'A broad folded grid with slow mirrored drift.',
    presets: [{ id: 'folded-grid', name: 'Folded Grid', values: {
      'pattern-freq': 3.565, speed: 0.235, 'pattern-mix': 1, drift: 1,
      'pole-fade': 1.4, 'camera-wander': 1, 'outer-rotation': 0.29530972,
      'outer-cell-x': 5.381125, 'outer-offset-x': 1.344, 'outer-offset-y': -1.456,
      'edge-width': 0.5, 'hue-shift-amount': 0.424,
      'hue-noise-scale': 2.2033439, 'palette-chroma': 0.4,
    } }],
  },
  {
    id: 'glitch-grid', display: 'Glitch Grid', source: 'grid', sourceKey: 'grid',
    projection: 'gnomonic-folded', lens: 'glitch', outer: 'mirror-tile', outerKey: 'mirror',
    inner: 'identity', innerKey: 'none', transfer: 'linear', coverage: 'edge-fade', valueKey: 'edge',
    palette: 'triadic', hue: 'noise', brightness: 'none', animatedProjection: false,
    description: 'A mirrored grid folded by the glitch lens.',
    presets: [{ id: 'folded-glitch', name: 'Folded Glitch', values: {
      'pattern-freq': 3.565, speed: 0.235, 'pattern-mix': 1, drift: 1,
      'pole-fade': 1.4, 'camera-wander': 1, 'outer-rotation': 0.29530972,
      'outer-cell-x': 5.381125, 'outer-offset-x': 1.344, 'outer-offset-y': -1.456,
      'edge-width': 0.5,
    } }],
  },
  {
    id: 'facet-wave', display: 'Facet Wave', source: 'grid', sourceKey: 'grid',
    projection: 'gnomonic-folded', lens: 'dodecahedral-kaleidoscope',
    outer: 'wave-shear', outerKey: 'wave', inner: 'mirror-tile', innerKey: 'mirror',
    transfer: 'linear', coverage: 'projection-squared', palette: 'triadic', hue: 'noise',
    brightness: 'none', animatedProjection: false,
    description: 'A wave-sheared grid moving across dodecahedral facets.',
    presets: [
      { id: 'wave-mirror', name: 'Wave Mirror', values: {
        'pattern-freq': 6.3287, speed: 0.04, complexity: 1.704, drift: 0.8,
        'source-angle-speed': 0.027, 'pole-fade': 2.311, 'camera-wander': 1,
        'outer-strength': -0.176, 'outer-speed': -0.00325,
        'outer-frequency': 1.408, 'outer-field-angle': 2.2305307,
        'hue-shift-amount': 0.721, 'palette-chroma': 1,
      } },
      { id: 'cup-hue', name: 'Cup Hue', values: {
        'pattern-freq': 6.3287, speed: 0.04, complexity: 1.704, drift: 0.8,
        'source-angle-speed': 0.027, 'pole-fade': 2.311, 'camera-wander': 1,
        'outer-strength': -0.176, 'outer-speed': -0.00325,
        'outer-frequency': 1.408, 'outer-field-angle': 2.2305307,
        'hue-shift-amount': 1, 'hue-noise-scale': 1.9717969,
        'palette-chroma': 1, 'palette-mapping': 'cup',
      } },
    ],
  },
  {
    id: 'contour-lattice', display: 'Contour Lattice', source: 'primitive-lattice', sourceKey: 'lattice',
    projection: 'gnomonic-folded', lens: 'identity', outer: 'affine-frame', outerKey: 'affine',
    inner: 'identity', innerKey: 'none', transfer: 'iso-contour', coverage: 'projection', valueKey: 'iso',
    palette: 'triadic', hue: 'noise', brightness: 'none', animatedProjection: true,
    description: 'An affine primitive lattice rendered as soft contours.',
    presets: [{ id: 'affine-contour', name: 'Affine Contour', values: {
      'lattice-cell-scale': 1.22925, 'lattice-shape': 1,
      'lattice-softness': 0.1608203, 'lattice-radius': 0.332981884,
      'projection-spin-speed': 0.0208791979, 'projection-wander': 0.00309175253,
      'camera-wander': 1, 'outer-speed': 0.015625,
      'outer-translation-x': 4, 'outer-translation-y': 4,
      'iso-level': 0.138, 'iso-width': 0.227034181,
      'hue-shift-amount': 0.398, 'hue-noise-scale': 0.8300313,
      'hue-noise-speed': 0.000212000014,
    } }],
  },
  {
    id: 'prism-lattice', display: 'Prism Lattice', source: 'primitive-lattice', sourceKey: 'lattice',
    projection: 'stereographic', lens: 'pentagonal-prism-kaleidoscope',
    outer: 'polar-chart-linear', outerKey: 'polar', inner: 'wave-shear', innerKey: 'wave',
    transfer: 'linear', coverage: 'projection-squared', palette: 'analogous', hue: 'noise',
    brightness: 'none', animatedProjection: true,
    description: 'A polar lattice folded through a pentagonal prism.',
    presets: [{ id: 'polar-wave', name: 'Polar Wave', values: {
      'lattice-cell-scale': 0.774140596, 'lattice-shape': 1,
      'lattice-softness': 0.377608389, 'lattice-radius': 0.290762514,
      'pole-fade': 2.273, 'projection-wander': 1, 'camera-wander': 1,
      'outer-speed': 0.000343749998, 'inner-speed': 0.000999999931,
      'hue-shift-amount': 0.268000007, 'hue-noise-scale': 2,
      'palette-chroma': 1, 'mapping-phase': -0.165999994,
      'palette-mapping': 'cup',
    } }],
  },
  {
    id: 'vector-facets', display: 'Vector Facets', source: 'grid', sourceKey: 'grid',
    projection: 'gnomonic-folded', lens: 'dodecahedral-kaleidoscope',
    outer: 'vector-noise-simplex', outerKey: 'vector', inner: 'mirror-tile', innerKey: 'mirror',
    transfer: 'linear', coverage: 'projection-squared', palette: 'triadic', hue: 'noise',
    brightness: 'cup', animatedProjection: false,
    description: 'A vector-noise grid refracted across dodecahedral facets.',
    presets: [{ id: 'vector-mirror', name: 'Vector Mirror', values: {
      'pattern-freq': 4.9755, speed: 0.04, complexity: 1.704, drift: 0.8,
      'source-angle-speed': 0.027, 'pole-fade': 2.311, 'camera-wander': 1,
      'outer-strength': 0.138, 'outer-speed': -0.00005,
      'inner-speed': 0.00327999983, 'hue-shift-amount': 0.721,
      'palette-chroma': 1, 'brightness-depth': 0.655, 'palette-mapping': 'cup',
    } }],
  },
  {
    id: 'hex-wave', display: 'Hex Wave', source: 'twin-wave', sourceKey: 'twin',
    projection: 'stereographic', lens: 'hexagonal-prism-kaleidoscope',
    outer: 'identity', outerKey: 'none', inner: 'mirror-tile', innerKey: 'mirror',
    transfer: 'linear', coverage: 'projection-squared', palette: 'analogous', hue: 'noise',
    brightness: 'none', animatedProjection: true,
    description: 'A twin wave folded through a hexagonal prism.',
    presets: [
      { id: 'hex-twin-wave', name: 'Hex Twin Wave', values: {
        'pattern-freq': 3.881, speed: 0.128598228, drift: 0.8,
        'source-angle-speed': 0.027, 'pole-fade': 4.971,
        'projection-wander': 1, 'camera-wander': 1,
        'hue-shift-amount': 0.226, 'hue-noise-scale': 1.47215629,
        'hue-noise-speed': 0.000138, 'palette-chroma': 1,
        'mapping-frequency': 1.341, 'mapping-phase': -1, 'palette-mapping': 'bell',
      } },
      { id: 'hex-twin-wave-alt', name: 'Hex Twin Wave Alt', values: {
        'pattern-freq': 3.881, speed: 0.12859823, drift: 0.8,
        'source-angle-speed': 0.027, 'pole-fade': 4.971,
        'projection-wander': 1, 'camera-wander': 1,
        'hue-shift-amount': 0.226, 'hue-noise-scale': 1.4721563,
        'hue-noise-speed': 0.000138, 'palette-chroma': 1,
        'mapping-frequency': 2, 'mapping-phase': -1, 'palette-mapping': 'bell',
      } },
    ],
  },
  {
    id: 'equator-grid', display: 'Equator Grid', source: 'grid', sourceKey: 'grid',
    projection: 'equirectangular', lens: 'dodecahedral-kaleidoscope',
    outer: 'identity', outerKey: 'none', inner: 'mirror-tile', innerKey: 'mirror',
    transfer: 'linear', coverage: 'projection-squared', palette: 'analogous', hue: 'noise',
    brightness: 'none', animatedProjection: true, centralMeridian: true,
    description: 'Dodecahedral grids mapped continuously around the equator.',
    presets: [
      { id: 'double-map', name: 'Double Map', values: {
        'pattern-freq': 3.9407, complexity: 3, 'pattern-mix': 1, drift: 0.8,
        'source-angle-speed': 0.0269999988, 'pole-fade': 2.14,
        'projection-wander': 0.165, 'camera-wander': 1,
        'inner-speed': 0.00013, 'inner-cell-y': 0.997703135,
        'hue-shift-amount': 0.366, 'hue-noise-scale': 1.47215629,
        'palette-chroma': 1, 'mapping-frequency': 2, 'palette-mapping': 'cup',
      } },
      { id: 'open-grid', name: 'Open Grid', values: {
        'pattern-freq': 3.9407, complexity: 3, 'pattern-mix': 1, drift: 0.8,
        'source-angle-speed': 0.0269999988, 'pole-fade': 2.14,
        'projection-wander': 0.165, 'camera-wander': 1,
        'inner-speed': 0.00013, 'inner-cell-y': 0.997703135,
        'hue-shift-amount': 0.366, 'hue-noise-scale': 1.47215629,
        'palette-chroma': 1, 'mapping-frequency': 1, 'palette-mapping': 'cup',
      } },
      { id: 'fine-grid', name: 'Fine Grid', values: {
        'pattern-freq': 0.3985, complexity: 3, 'pattern-mix': 1, drift: 0.8,
        'source-angle-speed': 0.0269999988, 'pole-fade': 2.14,
        'projection-wander': 0.165, 'camera-wander': 1,
        'inner-speed': 0.00058, 'inner-cell-y': 0.901890635,
        'hue-shift-amount': 0.366, 'hue-noise-scale': 1.47215629,
        'palette-chroma': 1, 'mapping-frequency': 21.212, 'palette-mapping': 'cup',
      } },
    ],
  },
  {
    id: 'cosmic-eyeball', display: 'Cosmic Eyeball', source: 'grid', sourceKey: 'grid',
    projection: 'stereographic', lens: 'glitch', outer: 'mirror-tile', outerKey: 'mirror',
    inner: 'identity', innerKey: 'none', transfer: 'linear', coverage: 'edge-fade', valueKey: 'edge',
    palette: 'triadic', hue: 'path-length', brightness: 'none', animatedProjection: false,
    description: 'A high-contrast mirrored grid with displacement-driven hue.',
    presets: [{ id: 'mirrored-grid', name: 'Mirrored Grid', values: {
      'pattern-freq': 2.5477, speed: 0.235, complexity: 1.854, drift: 1,
      'pole-fade': 1.4, 'camera-wander': 1, 'outer-rotation': 0.295309722,
      'outer-cell-x': 5.381125, 'outer-offset-x': 1.344, 'outer-offset-y': -1.456,
      'edge-width': 0.5, 'hue-shift-amount': 2.048, 'palette-chroma': 0.292,
    } }],
  },
  {
    id: 'mobius-grid', display: 'Mobius Grid', source: 'twin-wave', sourceKey: 'twin',
    projection: 'stereographic', lens: 'mobius', outer: 'identity', outerKey: 'none',
    inner: 'mirror-tile', innerKey: 'mirror', transfer: 'linear', coverage: 'projection-squared',
    palette: 'complementary', hue: 'path-length', brightness: 'cup', animatedProjection: true,
    description: 'A continuously animated Mobius lens over a mirrored twin wave.',
    presets: [
      { id: 'mobius-grid', name: 'Mobius Grid', values: {
        'pattern-freq': 10.158, speed: 0.245, drift: 0.8, 'source-angle-speed': 0.027,
        'pole-fade': 2.102, 'projection-wander': 1, 'camera-wander': 1,
        'mobius-a-re': -1.072, 'mobius-a-im': 0.304, 'mobius-b-re': 0.416,
        'mobius-b-im': 0, 'mobius-c-re': 0, 'mobius-c-im': 0,
        'mobius-d-re': 0.70710677, 'mobius-d-im': 0,
        'hue-shift-amount': 0.312, 'palette-chroma': 0.398,
      } },
      { id: 'mobius-grid-2', name: 'Mobius Grid 2', values: {
        'pattern-freq': 10.158, speed: 0.245, drift: 0.8, 'source-angle-speed': 0.027,
        'pole-fade': 2.102, 'projection-wander': 1, 'camera-wander': 1,
        'inner-speed': 0.005875, 'inner-cell-x': 0.2791094, 'inner-cell-y': 6.810328,
        'mobius-a-re': -1.072, 'mobius-a-im': 0.304, 'mobius-b-re': 0.416,
        'mobius-b-im': 0, 'mobius-c-re': 0, 'mobius-c-im': 0,
        'mobius-d-re': 0.70710677, 'mobius-d-im': 0,
        'hue-shift-amount': 0.312, 'palette-chroma': 0.398,
      } },
    ],
  },
];

// The specs above stay in the v1 six-role shape; compileShaderDocument runs
// them through expandV1Document against the operator catalog — the single code
// path every schema_version 1 document shares — so the committed output is the
// expansion's canonical v2 export.
const catalog = JSON.parse(
  await readFile(resolve(ROOT, 'scripts', 'engine_catalog.json'), 'utf8'));
const stale = [];
for (const spec of effects) {
  const compiled = compileShaderDocument(documentFor(spec), { catalog });
  if (compiled.status !== 'VALID') {
    console.error(`${spec.id}:`, JSON.stringify(compiled.diagnostics, null, 2));
    process.exit(1);
  }
  const name = `${spec.id.replaceAll('-', '_')}.shader.json`;
  const output = resolve(ROOT, 'patterns', name);
  const json = exportShaderDocumentJson(compiled.document);
  if (!CHECK) {
    await writeFile(output, json, 'utf8');
    continue;
  }
  // The committed blobs are LF, while core.autocrlf hands a Windows checkout a
  // CRLF working copy; compare the LF content, not the bytes on disk.
  const committed = await readFile(output, 'utf8').then(
    (text) => text.replaceAll('\r\n', '\n'), () => null);
  if (committed !== json) {
    stale.push(committed === null ? `${name} (missing)` : name);
  }
}

if (CHECK) {
  if (stale.length) {
    console.error('::error::patterns/ is out of sync with '
      + 'scripts/generate_promoted_shader_documents.mjs');
    for (const name of stale) {
      console.error(`  patterns/${name}`);
    }
    console.error(
      'Regenerate with: node scripts/generate_promoted_shader_documents.mjs');
    process.exit(1);
  }
  const noncanonical = [];
  const patternNames = (await readdir(resolve(ROOT, 'patterns')))
    .filter((name) => name.endsWith('.shader.json'));
  for (const name of patternNames) {
    const source = await readFile(resolve(ROOT, 'patterns', name), 'utf8')
      .then((text) => text.replaceAll('\r\n', '\n'));
    const compiled = compileShaderDocument(source, { catalog });
    if (compiled.status !== 'VALID' || exportShaderDocumentJson(compiled.document) !== source)
      noncanonical.push(name);
  }
  if (noncanonical.length) {
    console.error('::error::patterns/ contains noncanonical shader documents');
    for (const name of noncanonical)
      console.error(`  patterns/${name}`);
    process.exit(1);
  }
  console.log(
    `patterns/ matches the generator in full (${effects.length} documents) `
      + `and is canonical (${patternNames.length} documents).`);
}
