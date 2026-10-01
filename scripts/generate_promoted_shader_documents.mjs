import { exitAfterStderr } from './exit.mjs';
import { readFile, writeFile } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  compilePatternDocuments,
  loadOperatorCatalog,
} from './pattern_documents.mjs';
import {
  compileShaderDocument,
  declarationFromCatalogField,
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
  await exitAfterStderr(2);
}

// Choreography mirrors the shared ComposedEffect cadence.
const PRESET_DWELL_FRAMES = 600;
const INHERITED_SEGUE_FRAMES = 480;

const effects = [
  {
    id: "alien-brain",
    metadata: {
      "description": "Glitch-folded grids pulled through an animated wave shear.",
      "display_name": "Alien Brain",
    },
    chain: [
      {"label": "camera", "operator": "sphere.rotate.v2"},
      {"label": "lens", "operator": "sphere.lens.glitch.v2"},
      {"label": "project", "operator": "project.stereographic.v2"},
      {"label": "warp1", "operator": "warp.wave-shear.v2"},
      {"label": "sample", "operator": "sample.grid.v3"},
      {"label": "colorize", "operator": "colorize.generated-palette.v3"},
    ],
    values: {
      "camera.wander": 0,
      "colorize.brightness-envelope": "none",
      "colorize.hue-noise-scale": 1,
      "colorize.hue-noise-speed": 0,
      "colorize.hue-shift-amount": 0,
      "colorize.hue-shift-mode": "noise",
      "colorize.mapping-frequency": 1,
      "colorize.mapping-phase": 0,
      "colorize.palette-chroma": 0.6200000047683716,
      "colorize.palette-mapping": "linear",
      "colorize.palette-mode": "triadic",
      "colorize.phase-oscillation-depth": 0,
      "colorize.phase-oscillation-speed": 0,
      "colorize.value-opacity-high": 1,
      "colorize.value-opacity-low": 1,
      "project.projection-spin-speed": 0,
      "project.projection-wander": 0,
      "project.singularity-fade": 1,
      "sample.angle-speed": 0,
      "sample.complexity": 0,
      "sample.coverage-mode": "weight-squared",
      "sample.drift": 0,
      "sample.pattern-freq": 1,
      "sample.pattern-mix": 0,
      "sample.speed": 0,
      "sample.weight-mode": "projection",
      "warp1.field-angle": 0,
      "warp1.frequency": 1,
      "warp1.speed": 0,
      "warp1.strength": 0,
    },
    snap: [],
    presets: [
      {"id": "alien-brain", "name": "Alien Brain", "values": {"camera.wander": 0.800000011920929, "colorize.hue-noise-scale": 0.6304218769073486, "colorize.hue-shift-amount": 0.2919999957084656, "colorize.palette-chroma": 0.7879999876022339, "colorize.palette-mapping": "cup", "sample.complexity": 0.5, "sample.pattern-freq": 4.439000129699707, "sample.speed": 0.24500000476837158, "warp1.speed": 0.015625, "warp1.strength": 0.5}},
      {"id": "alien-brain-2", "name": "Alien Brain 2", "values": {"camera.wander": 0.800000011920929, "colorize.hue-noise-scale": 0.6304218769073486, "colorize.hue-shift-amount": 0.2919999957084656, "colorize.palette-chroma": 0.7879999876022339, "colorize.palette-mapping": "cup", "sample.complexity": 0.5, "sample.pattern-freq": 3.144700050354004, "sample.speed": 0.24500000476837158, "warp1.speed": 0.006906250026077032, "warp1.strength": 2.7200000286102295}},
      {"id": "alien-brain-3", "name": "Alien Brain 3", "values": {"camera.wander": 0.800000011920929, "colorize.hue-noise-scale": 0.6304218769073486, "colorize.hue-shift-amount": 0.2919999957084656, "colorize.palette-chroma": 0.7879999876022339, "colorize.palette-mapping": "cup", "sample.complexity": 1.6979999542236328, "sample.pattern-freq": 7.52269983291626, "sample.speed": 0.24500000476837158, "warp1.speed": 0.006906250026077032}},
      {"id": "alien-brain-4", "name": "Alien Brain 4", "values": {"camera.wander": 0.800000011920929, "colorize.hue-noise-scale": 0.6304218769073486, "colorize.hue-shift-amount": 0.2919999957084656, "colorize.palette-chroma": 0.7879999876022339, "colorize.palette-mapping": "cup", "sample.complexity": 1.6979999542236328, "sample.pattern-freq": 8.816200256347656, "sample.speed": 0.24500000476837158, "warp1.speed": 0.005593750160187483, "warp1.strength": 1.3760000467300415}},
    ],
  },
  {
    id: "kaleidoscope-hex-soft",
    metadata: {
      "description": "A drifting twin wave reflected through a kaleidoscope.",
      "display_name": "Kaleidoscope Hex Soft",
    },
    chain: [
      {"label": "camera", "operator": "sphere.rotate.v2"},
      {"label": "lens", "operator": "sphere.lens.kaleidoscope.v2"},
      {"label": "project", "operator": "project.stereographic.v2"},
      {"label": "warp2", "operator": "warp.mirror-tile.v2"},
      {"label": "sample", "operator": "sample.twin-wave.v3"},
      {"label": "colorize", "operator": "colorize.generated-palette.v3"},
    ],
    values: {
      "camera.wander": 0,
      "colorize.brightness-envelope": "none",
      "colorize.hue-noise-scale": 1,
      "colorize.hue-noise-speed": 0,
      "colorize.hue-shift-amount": 0,
      "colorize.hue-shift-mode": "noise",
      "colorize.mapping-frequency": 1,
      "colorize.mapping-phase": 0,
      "colorize.palette-chroma": 0.6200000047683716,
      "colorize.palette-mapping": "linear",
      "colorize.palette-mode": "triadic",
      "colorize.phase-oscillation-depth": 0,
      "colorize.phase-oscillation-speed": 0,
      "colorize.value-opacity-high": 1,
      "colorize.value-opacity-low": 1,
      "lens.symmetry": "azimuthal",
      "project.projection-spin-speed": 0,
      "project.projection-wander": 0,
      "project.singularity-fade": 1,
      "sample.angle-speed": 0,
      "sample.coverage-mode": "weight-squared",
      "sample.drift": 0,
      "sample.pattern-freq": 1,
      "sample.speed": 0,
      "sample.weight-mode": "projection",
      "warp2.cell-x": 1,
      "warp2.cell-y": 1,
      "warp2.offset-x": 0,
      "warp2.offset-y": 0,
      "warp2.rotation": 0,
      "warp2.speed": 0,
    },
    snap: [],
    presets: [
      {"id": "twin-wave", "name": "Twin Wave", "values": {"camera.wander": 1, "colorize.hue-noise-scale": 2.2033438682556152, "colorize.hue-noise-speed": -0.00040800002170726657, "colorize.hue-shift-amount": 0.27000001072883606, "colorize.palette-chroma": 0.3610000014305115, "project.projection-wander": 1, "project.singularity-fade": 4.9710001945495605, "sample.angle-speed": 0.05000000074505806, "sample.drift": 0.800000011920929, "sample.pattern-freq": 4.975500106811523, "sample.speed": 0.125}},
    ],
  },
  {
    id: "alien-ocean",
    metadata: {
      "description": "A broad folded grid with slow mirrored drift.",
      "display_name": "Alien Ocean",
    },
    chain: [
      {"label": "camera", "operator": "sphere.rotate.v2"},
      {"label": "lens", "operator": "sphere.lens.kaleidoscope.v2"},
      {"label": "project", "operator": "project.gnomonic.v2"},
      {"label": "warp1", "operator": "warp.mirror-tile.v2"},
      {"label": "sample", "operator": "sample.grid.v3"},
      {"label": "colorize", "operator": "colorize.generated-palette.v3"},
    ],
    values: {
      "camera.wander": 0,
      "colorize.brightness-envelope": "none",
      "colorize.hue-noise-scale": 1,
      "colorize.hue-noise-speed": 0,
      "colorize.hue-shift-amount": 0,
      "colorize.hue-shift-mode": "noise",
      "colorize.mapping-frequency": 1,
      "colorize.mapping-phase": 0,
      "colorize.palette-chroma": 0.6200000047683716,
      "colorize.palette-mapping": "linear",
      "colorize.palette-mode": "triadic",
      "colorize.phase-oscillation-depth": 0,
      "colorize.phase-oscillation-speed": 0,
      "colorize.value-opacity-high": 1,
      "colorize.value-opacity-low": 1,
      "lens.symmetry": "azimuthal",
      "project.frame": "identity",
      "project.hemisphere": "folded",
      "project.singularity-fade": 1,
      "sample.angle-speed": 0,
      "sample.complexity": 0,
      "sample.coverage-mode": "edge-fade",
      "sample.drift": 0,
      "sample.edge-width": 0.10000000149011612,
      "sample.pattern-freq": 1,
      "sample.pattern-mix": 0,
      "sample.speed": 0,
      "sample.weight-mode": "projection",
      "warp1.cell-x": 1,
      "warp1.cell-y": 1,
      "warp1.offset-x": 0,
      "warp1.offset-y": 0,
      "warp1.rotation": 0,
      "warp1.speed": 0,
    },
    snap: [],
    presets: [
      {"id": "folded-grid", "name": "Folded Grid", "values": {"camera.wander": 1, "colorize.hue-noise-scale": 2.2033438682556152, "colorize.hue-shift-amount": 0.42399999499320984, "colorize.palette-chroma": 0.4000000059604645, "project.singularity-fade": 1.399999976158142, "sample.drift": 1, "sample.edge-width": 0.5, "sample.pattern-freq": 3.565000057220459, "sample.pattern-mix": 1, "sample.speed": 0.23499999940395355, "warp1.cell-x": 5.381124973297119, "warp1.offset-x": 1.343999981880188, "warp1.offset-y": -1.4559999704360962, "warp1.rotation": 0.29530972242355347}},
    ],
  },
  {
    id: "alien-core",
    metadata: {
      "description": "A mirrored grid folded by the glitch lens.",
      "display_name": "Alien Core",
    },
    chain: [
      {"label": "camera", "operator": "sphere.rotate.v2"},
      {"label": "lens", "operator": "sphere.lens.glitch.v2"},
      {"label": "project", "operator": "project.gnomonic.v2"},
      {"label": "warp1", "operator": "warp.mirror-tile.v2"},
      {"label": "sample", "operator": "sample.grid.v3"},
      {"label": "colorize", "operator": "colorize.generated-palette.v3"},
    ],
    values: {
      "camera.wander": 0,
      "colorize.brightness-envelope": "none",
      "colorize.hue-noise-scale": 1,
      "colorize.hue-noise-speed": 0,
      "colorize.hue-shift-amount": 0,
      "colorize.hue-shift-mode": "noise",
      "colorize.mapping-frequency": 1,
      "colorize.mapping-phase": 0,
      "colorize.palette-chroma": 0.6200000047683716,
      "colorize.palette-mapping": "linear",
      "colorize.palette-mode": "triadic",
      "colorize.phase-oscillation-depth": 0,
      "colorize.phase-oscillation-speed": 0,
      "colorize.value-opacity-high": 1,
      "colorize.value-opacity-low": 1,
      "project.frame": "identity",
      "project.hemisphere": "folded",
      "project.singularity-fade": 1,
      "sample.angle-speed": 0,
      "sample.complexity": 0,
      "sample.coverage-mode": "edge-fade",
      "sample.drift": 0,
      "sample.edge-width": 0.10000000149011612,
      "sample.pattern-freq": 1,
      "sample.pattern-mix": 0,
      "sample.speed": 0,
      "sample.weight-mode": "projection",
      "warp1.cell-x": 1,
      "warp1.cell-y": 1,
      "warp1.offset-x": 0,
      "warp1.offset-y": 0,
      "warp1.rotation": 0,
      "warp1.speed": 0,
    },
    snap: [],
    presets: [
      {"id": "folded-glitch", "name": "Folded Glitch", "values": {"camera.wander": 1, "project.singularity-fade": 1.399999976158142, "sample.drift": 1, "sample.edge-width": 0.5, "sample.pattern-freq": 3.565000057220459, "sample.pattern-mix": 1, "sample.speed": 0.23499999940395355, "warp1.cell-x": 5.381124973297119, "warp1.offset-x": 1.343999981880188, "warp1.offset-y": -1.4559999704360962, "warp1.rotation": 0.29530972242355347}},
    ],
  },
  {
    id: "kaleidoscope-mandala",
    metadata: {
      "description": "A wave-sheared grid moving across dodecahedral facets.",
      "display_name": "Kaleidoscope Mandala",
    },
    chain: [
      {"label": "camera", "operator": "sphere.rotate.v2"},
      {"label": "lens", "operator": "sphere.lens.kaleidoscope.v2"},
      {"label": "project", "operator": "project.gnomonic.v2"},
      {"label": "warp1", "operator": "warp.wave-shear.v2"},
      {"label": "warp2", "operator": "warp.mirror-tile.v2"},
      {"label": "sample", "operator": "sample.grid.v3"},
      {"label": "colorize", "operator": "colorize.generated-palette.v3"},
    ],
    values: {
      "camera.wander": 0,
      "colorize.brightness-envelope": "none",
      "colorize.hue-noise-scale": 1,
      "colorize.hue-noise-speed": 0,
      "colorize.hue-shift-amount": 0,
      "colorize.hue-shift-mode": "noise",
      "colorize.mapping-frequency": 1,
      "colorize.mapping-phase": 0,
      "colorize.palette-chroma": 0.6200000047683716,
      "colorize.palette-mapping": "linear",
      "colorize.palette-mode": "triadic",
      "colorize.phase-oscillation-depth": 0,
      "colorize.phase-oscillation-speed": 0,
      "colorize.value-opacity-high": 1,
      "colorize.value-opacity-low": 1,
      "lens.symmetry": "dodecahedral",
      "project.frame": "identity",
      "project.hemisphere": "folded",
      "project.singularity-fade": 1,
      "sample.angle-speed": 0,
      "sample.complexity": 0,
      "sample.coverage-mode": "weight-squared",
      "sample.drift": 0,
      "sample.pattern-freq": 1,
      "sample.pattern-mix": 0,
      "sample.speed": 0,
      "sample.weight-mode": "projection",
      "warp1.field-angle": 0,
      "warp1.frequency": 1,
      "warp1.speed": 0,
      "warp1.strength": 0,
      "warp2.cell-x": 1,
      "warp2.cell-y": 1,
      "warp2.offset-x": 0,
      "warp2.offset-y": 0,
      "warp2.rotation": 0,
      "warp2.speed": 0,
    },
    snap: [],
    presets: [
      {"id": "wave-mirror", "name": "Wave Mirror", "values": {"camera.wander": 1, "colorize.hue-shift-amount": 0.7210000157356262, "colorize.palette-chroma": 1, "project.singularity-fade": 2.311000108718872, "sample.angle-speed": 0.027000000700354576, "sample.complexity": 1.7039999961853027, "sample.drift": 0.800000011920929, "sample.pattern-freq": 6.328700065612793, "sample.speed": 0.03999999910593033, "warp1.field-angle": 2.2305307388305664, "warp1.frequency": 1.4079999923706055, "warp1.speed": -0.0032500000670552254, "warp1.strength": -0.17599999904632568}},
      {"id": "cup-hue", "name": "Cup Hue", "values": {"camera.wander": 1, "colorize.hue-noise-scale": 1.9717968702316284, "colorize.hue-shift-amount": 1, "colorize.palette-chroma": 1, "colorize.palette-mapping": "cup", "project.singularity-fade": 2.311000108718872, "sample.angle-speed": 0.027000000700354576, "sample.complexity": 1.7039999961853027, "sample.drift": 0.800000011920929, "sample.pattern-freq": 6.328700065612793, "sample.speed": 0.03999999910593033, "warp1.field-angle": 2.2305307388305664, "warp1.frequency": 1.4079999923706055, "warp1.speed": -0.0032500000670552254, "warp1.strength": -0.17599999904632568}},
    ],
  },
  {
    id: "grid-space",
    metadata: {
      "description": "An affine primitive lattice rendered as soft contours.",
      "display_name": "Grid Space",
    },
    chain: [
      {"label": "camera", "operator": "sphere.rotate.v2"},
      {"label": "project", "operator": "project.gnomonic.v2"},
      {"label": "warp1", "operator": "warp.affine.v3"},
      {"label": "sample", "operator": "sample.lattice.v2"},
      {"label": "transfer", "operator": "field.transfer.iso-contour.v2"},
      {"label": "colorize", "operator": "colorize.generated-palette.v3"},
    ],
    values: {
      "camera.wander": 0,
      "colorize.brightness-envelope": "none",
      "colorize.hue-noise-scale": 1,
      "colorize.hue-noise-speed": 0,
      "colorize.hue-shift-amount": 0,
      "colorize.hue-shift-mode": "noise",
      "colorize.mapping-frequency": 1,
      "colorize.mapping-phase": 0,
      "colorize.palette-chroma": 0.6200000047683716,
      "colorize.palette-mapping": "linear",
      "colorize.palette-mode": "triadic",
      "colorize.phase-oscillation-depth": 0,
      "colorize.phase-oscillation-speed": 0,
      "colorize.value-opacity-high": 1,
      "colorize.value-opacity-low": 1,
      "project.hemisphere": "folded",
      "project.projection-spin-speed": 0,
      "project.projection-wander": 0,
      "project.singularity-fade": 1,
      "sample.coverage-mode": "weight",
      "sample.lattice-cell-scale": 1,
      "sample.lattice-radius": 0.25,
      "sample.lattice-shape": 0,
      "sample.lattice-softness": 0.05000000074505806,
      "sample.weight-mode": "projection",
      "transfer.iso-level": 0.5,
      "transfer.iso-width": 0.05000000074505806,
      "warp1.lattice-period": 1,
      "warp1.rotation-rate": 0,
      "warp1.scale-x": 1,
      "warp1.scale-y": 1,
      "warp1.shear": 0,
      "warp1.speed": 0,
      "warp1.translation-x": 0,
      "warp1.translation-y": 0,
    },
    snap: [],
    presets: [
      {"id": "affine-contour", "name": "Affine Contour", "values": {"camera.wander": 1, "colorize.hue-noise-scale": 0.8300312757492065, "colorize.hue-noise-speed": 0.00021200001356191933, "colorize.hue-shift-amount": 0.39800000190734863, "project.projection-spin-speed": 0.020879197865724564, "project.projection-wander": 0.003091752529144287, "sample.lattice-cell-scale": 1.2292499542236328, "sample.lattice-radius": 0.3329818844795227, "sample.lattice-shape": 1, "sample.lattice-softness": 0.16082030534744263, "transfer.iso-level": 0.1379999965429306, "transfer.iso-width": 0.22703418135643005, "warp1.lattice-period": 0.813504159450531, "warp1.speed": 0.015625, "warp1.translation-x": 4, "warp1.translation-y": 4}},
    ],
  },
  {
    id: "kaleidoscope-pent-bright",
    metadata: {
      "description": "A polar lattice folded through a pentagonal prism.",
      "display_name": "Kaleidoscope Pent Bright",
    },
    chain: [
      {"label": "camera", "operator": "sphere.rotate.v2"},
      {"label": "lens", "operator": "sphere.lens.kaleidoscope.v2"},
      {"label": "project", "operator": "project.stereographic.v2"},
      {"label": "warp1", "operator": "warp.polar-chart.v2"},
      {"label": "warp2", "operator": "warp.wave-shear.v2"},
      {"label": "sample", "operator": "sample.lattice.v2"},
      {"label": "colorize", "operator": "colorize.generated-palette.v3"},
    ],
    values: {
      "camera.wander": 0,
      "colorize.brightness-envelope": "none",
      "colorize.hue-noise-scale": 1,
      "colorize.hue-noise-speed": 0,
      "colorize.hue-shift-amount": 0,
      "colorize.hue-shift-mode": "noise",
      "colorize.mapping-frequency": 1,
      "colorize.mapping-phase": 0,
      "colorize.palette-chroma": 0.6200000047683716,
      "colorize.palette-mapping": "linear",
      "colorize.palette-mode": "analogous",
      "colorize.phase-oscillation-depth": 0,
      "colorize.phase-oscillation-speed": 0,
      "colorize.value-opacity-high": 1,
      "colorize.value-opacity-low": 1,
      "lens.symmetry": "pentagonal-prism",
      "project.projection-spin-speed": 0,
      "project.projection-wander": 0,
      "project.singularity-fade": 1,
      "sample.coverage-mode": "weight-squared",
      "sample.lattice-cell-scale": 1,
      "sample.lattice-radius": 0.25,
      "sample.lattice-shape": 0,
      "sample.lattice-softness": 0.05000000074505806,
      "sample.weight-mode": "projection",
      "warp1.angular-phase": 0,
      "warp1.mode": "linear",
      "warp1.radial-phase": 0,
      "warp1.radial-scale": 1,
      "warp1.speed": 0,
      "warp2.field-angle": 0,
      "warp2.frequency": 1,
      "warp2.speed": 0,
      "warp2.strength": 0,
    },
    snap: [],
    presets: [
      {"id": "polar-wave", "name": "Polar Wave", "values": {"camera.wander": 1, "colorize.hue-noise-scale": 2, "colorize.hue-shift-amount": 0.2680000066757202, "colorize.mapping-phase": -0.16599999368190765, "colorize.palette-chroma": 1, "colorize.palette-mapping": "cup", "project.projection-wander": 1, "project.singularity-fade": 2.2730000019073486, "sample.lattice-cell-scale": 0.7957746982574463, "sample.lattice-radius": 0.2907625138759613, "sample.lattice-shape": 1, "sample.lattice-softness": 0.37760838866233826, "warp1.speed": 0.00034374999813735485, "warp2.speed": 0.0009999999310821295}},
    ],
  },
  {
    id: "kaleidoscope-stained-glass",
    metadata: {
      "description": "A vector-noise grid refracted across dodecahedral facets.",
      "display_name": "Kaleidoscope Stained Glass",
    },
    chain: [
      {"label": "camera", "operator": "sphere.rotate.v2"},
      {"label": "lens", "operator": "sphere.lens.kaleidoscope.v2"},
      {"label": "project", "operator": "project.gnomonic.v2"},
      {"label": "warp1", "operator": "warp.vector-noise.v2"},
      {"label": "warp2", "operator": "warp.mirror-tile.v2"},
      {"label": "sample", "operator": "sample.grid.v3"},
      {"label": "colorize", "operator": "colorize.generated-palette.v3"},
    ],
    values: {
      "camera.wander": 0,
      "colorize.brightness-bottom": 0,
      "colorize.brightness-envelope": "cup",
      "colorize.brightness-top": 1,
      "colorize.hue-noise-scale": 1,
      "colorize.hue-noise-speed": 0,
      "colorize.hue-shift-amount": 0,
      "colorize.hue-shift-mode": "noise",
      "colorize.mapping-frequency": 1,
      "colorize.mapping-phase": 0,
      "colorize.palette-chroma": 0.6200000047683716,
      "colorize.palette-mapping": "linear",
      "colorize.palette-mode": "triadic",
      "colorize.phase-oscillation-depth": 0,
      "colorize.phase-oscillation-speed": 0,
      "colorize.value-opacity-high": 1,
      "colorize.value-opacity-low": 1,
      "lens.symmetry": "dodecahedral",
      "project.frame": "identity",
      "project.hemisphere": "folded",
      "project.singularity-fade": 1,
      "sample.angle-speed": 0,
      "sample.complexity": 0,
      "sample.coverage-mode": "weight-squared",
      "sample.drift": 0,
      "sample.pattern-freq": 1,
      "sample.pattern-mix": 0,
      "sample.speed": 0,
      "sample.weight-mode": "projection",
      "warp1.basis": "simplex",
      "warp1.scale": 1,
      "warp1.speed": 0,
      "warp1.strength": 0,
      "warp1.vector-angle": 0,
      "warp2.cell-x": 1,
      "warp2.cell-y": 1,
      "warp2.offset-x": 0,
      "warp2.offset-y": 0,
      "warp2.rotation": 0,
      "warp2.speed": 0,
    },
    snap: [],
    presets: [
      {"id": "vector-mirror", "name": "Vector Mirror", "values": {"camera.wander": 1, "colorize.brightness-bottom": 0.3449999988079071, "colorize.hue-shift-amount": 0.7210000157356262, "colorize.palette-chroma": 1, "colorize.palette-mapping": "cup", "project.singularity-fade": 2.311000108718872, "sample.angle-speed": 0.027000000700354576, "sample.complexity": 1.7039999961853027, "sample.drift": 0.800000011920929, "sample.pattern-freq": 4.975500106811523, "sample.speed": 0.03999999910593033, "warp1.speed": -4.999999873689376e-05, "warp1.strength": 0.1379999965429306, "warp2.speed": 0.003279999829828739}},
    ],
  },
  {
    id: "kaleidoscope-hex-bright",
    metadata: {
      "description": "A twin wave folded through a hexagonal prism.",
      "display_name": "Kaleidoscope Hex Bright",
    },
    chain: [
      {"label": "camera", "operator": "sphere.rotate.v2"},
      {"label": "lens", "operator": "sphere.lens.kaleidoscope.v2"},
      {"label": "project", "operator": "project.stereographic.v2"},
      {"label": "warp2", "operator": "warp.mirror-tile.v2"},
      {"label": "sample", "operator": "sample.twin-wave.v3"},
      {"label": "colorize", "operator": "colorize.generated-palette.v3"},
    ],
    values: {
      "camera.wander": 0,
      "colorize.brightness-envelope": "none",
      "colorize.hue-noise-scale": 1,
      "colorize.hue-noise-speed": 0,
      "colorize.hue-shift-amount": 0,
      "colorize.hue-shift-mode": "noise",
      "colorize.mapping-frequency": 1,
      "colorize.mapping-phase": 0,
      "colorize.palette-chroma": 0.6200000047683716,
      "colorize.palette-mapping": "linear",
      "colorize.palette-mode": "analogous",
      "colorize.phase-oscillation-depth": 0,
      "colorize.phase-oscillation-speed": 0,
      "colorize.value-opacity-high": 1,
      "colorize.value-opacity-low": 1,
      "lens.symmetry": "hexagonal-prism",
      "project.projection-spin-speed": 0,
      "project.projection-wander": 0,
      "project.singularity-fade": 1,
      "sample.angle-speed": 0,
      "sample.coverage-mode": "weight-squared",
      "sample.drift": 0,
      "sample.pattern-freq": 1,
      "sample.speed": 0,
      "sample.weight-mode": "projection",
      "warp2.cell-x": 1,
      "warp2.cell-y": 1,
      "warp2.offset-x": 0,
      "warp2.offset-y": 0,
      "warp2.rotation": 0,
      "warp2.speed": 0,
    },
    snap: [],
    presets: [
      {"id": "hex-twin-wave", "name": "Hex Twin Wave", "values": {"camera.wander": 1, "colorize.hue-noise-scale": 1.472156286239624, "colorize.hue-noise-speed": 0.00013800000306218863, "colorize.hue-shift-amount": 0.22599999606609344, "colorize.mapping-frequency": 1.340999960899353, "colorize.mapping-phase": -1, "colorize.palette-chroma": 1, "colorize.palette-mapping": "bell", "project.projection-wander": 1, "project.singularity-fade": 4.9710001945495605, "sample.angle-speed": 0.027000000700354576, "sample.drift": 0.800000011920929, "sample.pattern-freq": 3.88100004196167, "sample.speed": 0.12859822809696198}},
      {"id": "hex-twin-wave-alt", "name": "Hex Twin Wave Alt", "values": {"camera.wander": 1, "colorize.hue-noise-scale": 1.472156286239624, "colorize.hue-noise-speed": 0.00013800000306218863, "colorize.hue-shift-amount": 0.22599999606609344, "colorize.mapping-frequency": 2, "colorize.mapping-phase": -1, "colorize.palette-chroma": 1, "colorize.palette-mapping": "bell", "project.projection-wander": 1, "project.singularity-fade": 4.9710001945495605, "sample.angle-speed": 0.027000000700354576, "sample.drift": 0.800000011920929, "sample.pattern-freq": 3.88100004196167, "sample.speed": 0.12859822809696198}},
    ],
  },
  {
    id: "kaleidoscope-flowers",
    metadata: {
      "description": "Dodecahedral grids mapped continuously around the equator.",
      "display_name": "Kaleidoscope Flowers",
    },
    chain: [
      {"label": "camera", "operator": "sphere.rotate.v2"},
      {"label": "lens", "operator": "sphere.lens.kaleidoscope.v2"},
      {"label": "project", "operator": "project.equirectangular.v2"},
      {"label": "warp2", "operator": "warp.mirror-tile.v2"},
      {"label": "sample", "operator": "sample.grid.v3"},
      {"label": "colorize", "operator": "colorize.generated-palette.v3"},
    ],
    values: {
      "camera.wander": 0,
      "colorize.brightness-envelope": "none",
      "colorize.hue-noise-scale": 1,
      "colorize.hue-noise-speed": 0,
      "colorize.hue-shift-amount": 0,
      "colorize.hue-shift-mode": "noise",
      "colorize.mapping-frequency": 1,
      "colorize.mapping-phase": 0,
      "colorize.palette-chroma": 0.6200000047683716,
      "colorize.palette-mapping": "linear",
      "colorize.palette-mode": "analogous",
      "colorize.phase-oscillation-depth": 0,
      "colorize.phase-oscillation-speed": 0,
      "colorize.value-opacity-high": 1,
      "colorize.value-opacity-low": 1,
      "lens.symmetry": "dodecahedral",
      "project.central-meridian": 0,
      "project.projection-spin-speed": 0,
      "project.projection-wander": 0,
      "project.singularity-fade": 1,
      "sample.angle-speed": 0,
      "sample.complexity": 0,
      "sample.coverage-mode": "weight-squared",
      "sample.drift": 0,
      "sample.pattern-freq": 1,
      "sample.pattern-mix": 0,
      "sample.speed": 0,
      "sample.weight-mode": "projection",
      "warp2.cell-x": 1,
      "warp2.cell-y": 1,
      "warp2.offset-x": 0,
      "warp2.offset-y": 0,
      "warp2.rotation": 0,
      "warp2.speed": 0,
    },
    snap: [],
    presets: [
      {"id": "double-map", "name": "Double Map", "values": {"camera.wander": 1, "colorize.hue-noise-scale": 1.472156286239624, "colorize.hue-shift-amount": 0.3659999966621399, "colorize.mapping-frequency": 2, "colorize.palette-chroma": 1, "colorize.palette-mapping": "cup", "project.projection-wander": 0.16500000655651093, "project.singularity-fade": 2.140000104904175, "sample.angle-speed": 0.026999998837709427, "sample.complexity": 3, "sample.drift": 0.800000011920929, "sample.pattern-freq": 3.940700054168701, "sample.pattern-mix": 1, "warp2.cell-x": 1.0471975803375244, "warp2.cell-y": 0.9977031350135803, "warp2.speed": 0.00013000000035390258}},
      {"id": "open-grid", "name": "Open Grid", "values": {"camera.wander": 1, "colorize.hue-noise-scale": 1.472156286239624, "colorize.hue-shift-amount": 0.3659999966621399, "colorize.palette-chroma": 1, "colorize.palette-mapping": "cup", "project.projection-wander": 0.16500000655651093, "project.singularity-fade": 2.140000104904175, "sample.angle-speed": 0.026999998837709427, "sample.complexity": 3, "sample.drift": 0.800000011920929, "sample.pattern-freq": 3.940700054168701, "sample.pattern-mix": 1, "warp2.cell-x": 1.0471975803375244, "warp2.cell-y": 0.9977031350135803, "warp2.speed": 0.00013000000035390258}},
      {"id": "fine-grid", "name": "Fine Grid", "values": {"camera.wander": 1, "colorize.hue-noise-scale": 1.472156286239624, "colorize.hue-shift-amount": 0.3659999966621399, "colorize.mapping-frequency": 21.211999893188477, "colorize.palette-chroma": 1, "colorize.palette-mapping": "cup", "project.projection-wander": 0.16500000655651093, "project.singularity-fade": 2.140000104904175, "sample.angle-speed": 0.026999998837709427, "sample.complexity": 3, "sample.drift": 0.800000011920929, "sample.pattern-freq": 0.398499995470047, "sample.pattern-mix": 1, "warp2.cell-x": 1.0471975803375244, "warp2.cell-y": 0.9018906354904175, "warp2.speed": 0.0005799999926239252}},
    ],
  },
  {
    id: "cosmic-eyeball",
    metadata: {
      "description": "A high-contrast mirrored grid with displacement-driven hue.",
      "display_name": "Cosmic Eyeball",
    },
    chain: [
      {"label": "camera", "operator": "sphere.rotate.v2"},
      {"label": "lens", "operator": "sphere.lens.glitch.v2"},
      {"label": "project", "operator": "project.stereographic.v2"},
      {"label": "warp1", "operator": "warp.mirror-tile.v2"},
      {"label": "sample", "operator": "sample.grid.v3"},
      {"label": "colorize", "operator": "colorize.generated-palette.v3"},
    ],
    values: {
      "camera.wander": 0,
      "colorize.brightness-envelope": "none",
      "colorize.hue-shift-amount": 0,
      "colorize.hue-shift-mode": "path-length",
      "colorize.mapping-frequency": 1,
      "colorize.mapping-phase": 0,
      "colorize.palette-chroma": 0.6200000047683716,
      "colorize.palette-mapping": "linear",
      "colorize.palette-mode": "triadic",
      "colorize.phase-oscillation-depth": 0,
      "colorize.phase-oscillation-speed": 0,
      "colorize.value-opacity-high": 1,
      "colorize.value-opacity-low": 1,
      "project.frame": "identity",
      "project.singularity-fade": 1,
      "sample.angle-speed": 0,
      "sample.complexity": 0,
      "sample.coverage-mode": "edge-fade",
      "sample.drift": 0,
      "sample.edge-width": 0.10000000149011612,
      "sample.pattern-freq": 1,
      "sample.pattern-mix": 0,
      "sample.speed": 0,
      "sample.weight-mode": "projection",
      "warp1.cell-x": 1,
      "warp1.cell-y": 1,
      "warp1.offset-x": 0,
      "warp1.offset-y": 0,
      "warp1.rotation": 0,
      "warp1.speed": 0,
    },
    snap: [],
    presets: [
      {"id": "mirrored-grid", "name": "Mirrored Grid", "values": {"camera.wander": 1, "colorize.hue-shift-amount": 2.0480000972747803, "colorize.palette-chroma": 0.2919999957084656, "project.singularity-fade": 1.399999976158142, "sample.complexity": 1.8539999723434448, "sample.drift": 1, "sample.edge-width": 0.5, "sample.pattern-freq": 2.5476999282836914, "sample.speed": 0.23499999940395355, "warp1.cell-x": 5.381124973297119, "warp1.offset-x": 1.343999981880188, "warp1.offset-y": -1.4559999704360962, "warp1.rotation": 0.29530972242355347}},
    ],
  },
  {
    id: "mobius-grid",
    metadata: {
      "description": "A continuously animated Mobius lens over a mirrored twin wave.",
      "display_name": "Mobius Grid",
    },
    chain: [
      {"label": "camera", "operator": "sphere.rotate.v2"},
      {"label": "lens", "operator": "sphere.lens.mobius.v2"},
      {"label": "project", "operator": "project.stereographic.v2"},
      {"label": "warp2", "operator": "warp.mirror-tile.v2"},
      {"label": "sample", "operator": "sample.twin-wave.v3"},
      {"label": "colorize", "operator": "colorize.generated-palette.v3"},
    ],
    values: {
      "camera.wander": 0,
      "colorize.brightness-bottom": 0,
      "colorize.brightness-envelope": "cup",
      "colorize.brightness-top": 1,
      "colorize.hue-shift-amount": 0,
      "colorize.hue-shift-mode": "path-length",
      "colorize.mapping-frequency": 1,
      "colorize.mapping-phase": 0,
      "colorize.palette-chroma": 0.6200000047683716,
      "colorize.palette-mapping": "linear",
      "colorize.palette-mode": "complementary",
      "colorize.phase-oscillation-depth": 0,
      "colorize.phase-oscillation-speed": 0,
      "colorize.value-opacity-high": 1,
      "colorize.value-opacity-low": 1,
      "lens.mobius-a-im": 0,
      "lens.mobius-a-re": 0.7071067690849304,
      "lens.mobius-b-im": 0,
      "lens.mobius-b-re": 0,
      "lens.mobius-c-im": 0,
      "lens.mobius-c-re": 0,
      "lens.mobius-d-im": 0,
      "lens.mobius-d-re": 0.7071067690849304,
      "project.projection-spin-speed": 0,
      "project.projection-wander": 0,
      "project.singularity-fade": 1,
      "sample.angle-speed": 0,
      "sample.coverage-mode": "weight-squared",
      "sample.drift": 0,
      "sample.pattern-freq": 1,
      "sample.speed": 0,
      "sample.weight-mode": "projection",
      "warp2.cell-x": 1,
      "warp2.cell-y": 1,
      "warp2.offset-x": 0,
      "warp2.offset-y": 0,
      "warp2.rotation": 0,
      "warp2.speed": 0,
    },
    snap: ["lens.mobius-a-im", "lens.mobius-a-re", "lens.mobius-b-im", "lens.mobius-b-re", "lens.mobius-c-im", "lens.mobius-c-re", "lens.mobius-d-im", "lens.mobius-d-re"],
    presets: [
      {"id": "mobius-grid", "name": "Mobius Grid", "values": {"camera.wander": 1, "colorize.hue-shift-amount": 0.31200000643730164, "colorize.palette-chroma": 0.39800000190734863, "lens.mobius-a-im": 0.30399999022483826, "lens.mobius-a-re": -1.0720000267028809, "lens.mobius-b-re": 0.41600000858306885, "project.projection-wander": 1, "project.singularity-fade": 2.1019999980926514, "sample.angle-speed": 0.027000000700354576, "sample.drift": 0.800000011920929, "sample.pattern-freq": 10.157999992370605, "sample.speed": 0.24500000476837158}},
      {"id": "mobius-grid-2", "name": "Mobius Grid 2", "values": {"camera.wander": 1, "colorize.hue-shift-amount": 0.31200000643730164, "colorize.palette-chroma": 0.39800000190734863, "lens.mobius-a-im": 0.30399999022483826, "lens.mobius-a-re": -1.0720000267028809, "lens.mobius-b-re": 0.41600000858306885, "project.projection-wander": 1, "project.singularity-fade": 2.1019999980926514, "sample.angle-speed": 0.027000000700354576, "sample.drift": 0.800000011920929, "sample.pattern-freq": 10.157999992370605, "sample.speed": 0.24500000476837158, "warp2.cell-x": 0.279109388589859, "warp2.cell-y": 6.810328006744385, "warp2.speed": 0.005874999798834324}},
    ],
  },
];

const documentFor = (spec, catalog) => {
  const operators = new Map(catalog.operators.map((operator) => [operator.id, operator]));
  const slots = new Map(spec.chain.map((slot) => [slot.label, operators.get(slot.operator)]));
  const parameters = Object.entries(spec.values).map(([id, value]) => {
    const [label, fieldId] = id.split('.');
    const operator = slots.get(label);
    const field = operator?.params.find((entry) => entry.id === fieldId);
    if (!field) throw new Error(`Unknown generated parameter ${spec.id}/${id}`);
    return {
      ...declarationFromCatalogField(label, field, operator.id),
      default: value,
      ...(spec.snap.includes(id) ? { interpolation: { kind: 'SNAP' } } : {}),
    };
  });
  const presets = spec.presets.map((preset) => ({
    preset_id: preset.id,
    display_name: preset.name,
    values: { ...spec.values, ...preset.values },
  }));
  return {
    schema_version: 2,
    catalog_version: catalog.catalog_version,
    document_id: spec.id,
    effect_id: spec.id,
    descriptor: {
      chain: spec.chain,
      parameters,
      path_policies: [{ id: 'parallel', kind: 'PARALLEL' }],
      serialization: { schema_version: 1, fields: parameters.map((parameter) => parameter.id) },
    },
    preset_bank: {
      schema_version: 1,
      presets,
      edges: presets.length < 2 ? [] : presets.map((preset, index) => ({
        from: preset.preset_id,
        to: presets[(index + 1) % presets.length].preset_id,
        path_policy: 'parallel', easing: 'EASE_IN_OUT_SIN',
        duration: INHERITED_SEGUE_FRAMES,
      })),
      absent_edge_fallback: {
        manual: 'SNAP', automatic: 'REJECT', synchronized: 'SNAP',
        restore: 'SNAP', authoring: 'SNAP',
      },
      choreography: {
        generated_order: presets.map((preset) => preset.preset_id),
        dwell: Object.fromEntries(presets.map(
          (preset) => [preset.preset_id, PRESET_DWELL_FRAMES])),
      },
    },
    effect_metadata: spec.metadata,
  };
};

const catalog = await loadOperatorCatalog();
const stale = [];
for (const spec of effects) {
  const document = documentFor(spec, catalog);
  const compiled = compileShaderDocument(document, { catalog });
  if (compiled.status !== 'VALID') {
    console.error(`${spec.id}:`, JSON.stringify(compiled.diagnostics, null, 2));
    await exitAfterStderr(1);
  }
  const name = `${spec.id.replaceAll('-', '_')}.shader.json`;
  const output = resolve(ROOT, 'patterns', name);
  const json = exportShaderDocumentJson(compiled.document);
  if (!CHECK) {
    await writeFile(output, json, 'utf8');
    continue;
  }
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
    await exitAfterStderr(1);
  }
  const noncanonical = [];
  const patterns = await compilePatternDocuments(catalog);
  for (const { name, source, compiled } of patterns) {
    if (compiled.status !== 'VALID' || exportShaderDocumentJson(compiled.document) !== source)
      noncanonical.push(name);
  }
  if (noncanonical.length) {
    console.error('::error::patterns/ contains noncanonical shader documents');
    for (const name of noncanonical)
      console.error(`  patterns/${name}`);
    await exitAfterStderr(1);
  }
  console.log(
    `patterns/ generated subset matches its specs (${effects.length} documents) `
      + `and is canonical (${patterns.length} documents).`);
}
