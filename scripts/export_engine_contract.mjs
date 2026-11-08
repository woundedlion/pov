import { readFile, writeFile } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

const [modulePath] = process.argv.slice(2);
if (!modulePath) throw new Error('usage: export_engine_contract.mjs <module.js>');
const directory = dirname(resolve(modulePath));
const { default: createModule } = await import(pathToFileURL(resolve(modulePath)));
const module = await createModule();
const numericEnums = {};
const classMethods = {};
for (const [name, value] of Object.entries(module)) {
  if (typeof value !== 'function') continue;
  const members = Object.entries(value).filter(([, entry]) => entry instanceof value);
  if (members.length) {
    numericEnums[name] = Object.fromEntries(members.map(([key, entry]) => [key, entry.value]));
  } else if (value.prototype) {
    const methods = Object.keys(value.prototype).filter((key) => key !== 'constructor');
    if (methods.length) classMethods[name] = [...methods, 'delete'].sort();
  }
}
if (!(module.SHOW_FRAMES_PER_SECOND > 0) || !classMethods.HolosphereEngine)
  throw new Error('engine contract is incomplete');

/**
 * Shortest decimal that rounds to the same float32 as a bound engine float.
 * @param {string} name - Bound constant name.
 * @returns {number} The decimal value.
 */
const engineFloat = (name) => {
  const value = module[name];
  if (!Number.isFinite(value)) throw new Error(`engine contract is missing ${name}`);
  for (let digits = 1; digits <= 9; digits++) {
    const decimal = Number(value.toPrecision(digits));
    if (Math.fround(decimal) === value) return decimal;
  }
  return value;
};
const unswept = module.MeshOps?.getUnsweptOps?.();
if (!Array.isArray(unswept) || !unswept.every((op) => typeof op === 'string'))
  throw new Error('engine contract is missing MeshOps.getUnsweptOps');
const morphSweep = {
  bands: {
    truncate: {
      t: {
        min: engineFloat('MORPH_TRUNCATE_MIN'),
        max: engineFloat('MORPH_TRUNCATE_MAX'),
        excluded: {
          lower: engineFloat('MORPH_TRUNCATE_PINCH_LOWER'),
          upper: engineFloat('MORPH_TRUNCATE_PINCH_UPPER'),
        },
      },
    },
    chamfer: { t: { min: engineFloat('MORPH_CHAMFER_MIN'), max: engineFloat('MORPH_CHAMFER_MAX') } },
  },
  unswept: [...unswept].sort(),
};
const engine = new module.HolosphereEngine();
let effectIdentities;
try { effectIdentities = engine.getEffectIdentities(); }
finally { engine.delete?.(); }
if (!Array.isArray(effectIdentities) || effectIdentities.length === 0)
  throw new Error('engine identity contract is incomplete');
await writeFile(resolve(directory, 'engine_contract.mjs'), `// @ts-check
// Generated from the installed engine's bound exports.
export const FPS = ${module.SHOW_FRAMES_PER_SECOND};
/** @type {Record<string, Record<string, number>>} */
export const NUMERIC_ENUMS = ${JSON.stringify(numericEnums, null, 2)};
/** @type {Record<string, string[]>} */
export const CLASS_METHODS = ${JSON.stringify(classMethods, null, 2)};
/** @type {Array<{id: string, display_name: string}>} */
export const EFFECT_IDENTITIES = ${JSON.stringify(effectIdentities, null, 2)};
/**
 * Recipe-morph sweep coverage of the authored primitive ops. A band entry is
 * the inclusive parameter range a leg covers, less any excluded interval
 * (above the lower bound through the upper bound); an unswept op never sweeps;
 * any other primitive sweeps at every parameter it accepts.
 * @type {{bands: Record<string, Record<string, {min: number, max: number, excluded?: {lower: number, upper: number}}>>, unswept: string[]}}
 */
export const MORPH_SWEEP = ${JSON.stringify(morphSweep, null, 2)};
`);

let declarations = await readFile(resolve(directory, 'holosphere_bindings.d.ts'), 'utf8');
const handles = [...declarations.matchAll(/export interface (\w+) extends ClassHandle/g)]
  .map((match) => match[1]);
for (const handle of handles) {
  declarations = declarations.replace(new RegExp(`^(?!\\s*new\\()[^\\n]*\\): ${handle};`, 'gm'),
    (line) => line.replace(`: ${handle};`, `: ${handle} | null;`));
}
await writeFile(resolve(directory, 'holosphere_bindings.d.ts'), declarations);
await writeFile(resolve(directory, 'holosphere_wasm.d.ts'), `// Generated binding aliases and authored runtime extensions.
import type { MainModule } from './holosphere_bindings.js';
import type { RuntimeState, ModuleOptions } from './engine_types.js';
export * from './holosphere_bindings.js';
export * from './engine_types.js';
export type HolosphereModule = MainModule & RuntimeState;
export type MeshHandle = import('./holosphere_bindings.js').MeshOps;
export type MeshOpsStatics = MainModule['MeshOps'];
export default function createHolosphereModule(options?: ModuleOptions): Promise<HolosphereModule>;
`);
