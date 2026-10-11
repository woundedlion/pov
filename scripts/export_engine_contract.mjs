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
/** @type {Array<{id: string, display_name: string, legacy_ids: string[]}>} */
export const EFFECT_IDENTITIES = ${JSON.stringify(effectIdentities, null, 2)};
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
