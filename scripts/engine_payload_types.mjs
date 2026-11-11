// Checks the engine's JavaScript payloads against their TypeScript declarations.
//
// Usage: node scripts/engine_payload_types.mjs <path/to/holosphere_wasm.js>
//
// Embind emits holosphere_bindings.d.ts next to the module; each bound function
// returning a registered payload alias (`type <Name>Value = ...`) is called and
// its result validated against that alias, with targets/wasm/engine_types.d.ts
// resolving the authored payload schemas. The run fails when a payload carries
// an undeclared field, lacks a required one or has the wrong type, when a
// payload-returning function has no sample call, or when an engine_types.d.ts
// declaration is never reached by a payload.
import { readFileSync } from 'node:fs';
import { dirname, isAbsolute, join, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
import { exitAfterStderr, isMain } from './exit.mjs';

const ENGINE_TYPES = new URL('../targets/wasm/engine_types.d.ts', import.meta.url);

/** Loader extensions in engine_types.d.ts that describe no engine payload. */
export const LOADER_TYPES = new Set(['RuntimeState', 'ModuleOptions']);

/**
 * A shader chain whose snapshot carries a modulation recipe and whose field
 * operator reports a fold advisory.
 */
const SAMPLE_CHAIN = [
  { instance: 'field', operator: 'sphere.displace.potential.v1' },
  { instance: 'project', operator: 'project.stereographic.v2' },
  { instance: 'sample', operator: 'sample.grid.v4' },
  { instance: 'colorize', operator: 'colorize.generated-palette.v4' },
  { instance: 'modulation', operator: 'color.modulate.v1' },
];

/** Writes that give {@link SAMPLE_CHAIN}'s field operator a fold threshold. */
const SAMPLE_PARAMETERS = [
  { name: 'field.strength', value: 0.1 },
  { name: 'field.term-1-amplitude', value: 1 },
  { name: 'field.term-1-l', value: 6 },
  { name: 'field.term-1-m', value: 6 },
];

const PRIMITIVES = new Set(['number', 'string', 'boolean', 'object', 'unknown', 'null']);

const TOKEN = /\s*(?:(\.\.\.|=>|[{}[\]()<>|;:,?.=&])|("(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*')|(-?\d+(?:\.\d+)?)|([A-Za-z_$][\w$]*))/y;

/**
 * Splits TypeScript type text into tokens.
 * @param {string} text - Type source.
 * @returns {string[]} Punctuation, literal and identifier tokens.
 */
function tokenize(text) {
  const tokens = [];
  TOKEN.lastIndex = 0;
  while (TOKEN.lastIndex < text.length) {
    if (/^\s*$/.test(text.slice(TOKEN.lastIndex))) break;
    const start = TOKEN.lastIndex;
    const match = TOKEN.exec(text);
    if (!match) throw new Error(`cannot tokenize type at "${text.slice(start, start + 20)}"`);
    tokens.push(match[1] ?? match[2] ?? match[3] ?? match[4]);
  }
  return tokens;
}

/**
 * Parses TypeScript type text into a type node.
 * @param {string} text - A type expression, or object-type members when
 *   `members` is set.
 * @param {{members?: boolean}} [options] - Parse a member list instead.
 * @returns {object} The type node.
 * @throws {Error} On syntax outside the supported subset.
 */
export function parseType(text, { members = false } = {}) {
  const tokens = tokenize(text);
  let at = 0;
  const peek = () => tokens[at];
  const take = (expected) => {
    const token = tokens[at++];
    if (expected !== undefined && token !== expected)
      throw new Error(`expected "${expected}" but found "${token}" in "${text}"`);
    return token;
  };
  const list = (close, item) => {
    const items = [];
    while (peek() !== close) {
      items.push(item());
      if (peek() === ',') take(',');
    }
    take(close);
    return items;
  };
  const memberList = (close) => {
    const fields = [];
    while (peek() !== close) {
      if (peek() === 'readonly' && tokens[at + 1] !== ':' && tokens[at + 1] !== '?') take();
      const name = take();
      if (!/^[A-Za-z_$][\w$]*$/.test(name)) throw new Error(`bad member name "${name}" in "${text}"`);
      const optional = peek() === '?';
      if (optional) take('?');
      take(':');
      fields.push({ name, optional, type: union() });
      if (peek() === ';' || peek() === ',') take();
    }
    if (close !== undefined) take(close);
    return { kind: 'object', members: fields };
  };
  const primary = () => {
    const token = take();
    if (token === '(') {
      const inner = union();
      take(')');
      return inner;
    }
    if (token === '{') return memberList('}');
    if (token === '[') return { kind: 'tuple', elements: list(']', union) };
    if (/^["']/.test(token)) return { kind: 'literal', value: token.slice(1, -1) };
    if (/^-?\d/.test(token)) return { kind: 'literal', value: Number(token) };
    if (token === 'true' || token === 'false') return { kind: 'literal', value: token === 'true' };
    if (token === 'import') {
      take('(');
      take();
      take(')');
      take('.');
      return { kind: 'ref', name: take(), args: [] };
    }
    if (PRIMITIVES.has(token)) return { kind: 'primitive', name: token };
    if (!/^[A-Za-z_$]/.test(token)) throw new Error(`unexpected "${token}" in "${text}"`);
    let args = [];
    if (peek() === '<') {
      take('<');
      args = list('>', union);
    }
    return { kind: 'ref', name: token, args };
  };
  const postfix = () => {
    let node = primary();
    while (peek() === '[' && tokens[at + 1] === ']') {
      take('[');
      take(']');
      node = { kind: 'array', element: node };
    }
    return node;
  };
  const union = () => {
    if (peek() === '|') take('|');
    const options = [postfix()];
    while (peek() === '|') {
      take('|');
      options.push(postfix());
    }
    return options.length === 1 ? options[0] : { kind: 'union', options };
  };
  const node = members ? memberList(undefined) : union();
  if (at !== tokens.length) throw new Error(`trailing "${tokens.slice(at).join(' ')}" in "${text}"`);
  return node;
}

/**
 * Index of the bracket closing the one opened just before `start`.
 * @param {string} text - Source text.
 * @param {number} start - Index after the opening bracket.
 * @param {string} open - Opening bracket.
 * @param {string} close - Closing bracket.
 * @returns {number} Index of the matching close.
 */
function matching(text, start, open, close) {
  let depth = 1;
  for (let index = start; index < text.length; ++index) {
    if (text[index] === open) ++depth;
    else if (text[index] === close && --depth === 0) return index;
  }
  throw new Error(`unbalanced "${open}" at offset ${start}`);
}

/**
 * Top-level interface and type-alias declarations, parsed on first use.
 * @param {string} source - Declaration file text.
 * @returns {Map<string, {kind: 'interface'|'alias', extends: string[], body: string}>}
 *   Raw declarations by name.
 */
export function declarations(source) {
  const text = source.replace(/\/\*[\s\S]*?\*\//g, '').replace(/\/\/[^\n]*/g, '');
  const found = new Map();
  const head = /(?:^|\n)\s*(?:export\s+)?(interface|type)\s+(\w+)(<[^>{=]*>)?/g;
  let match;
  while ((match = head.exec(text))) {
    const [, keyword, name, generics] = match;
    if (generics) continue;
    let body;
    let parents = [];
    if (keyword === 'interface') {
      const open = text.indexOf('{', head.lastIndex);
      const between = text.slice(head.lastIndex, open).trim();
      if (between.startsWith('extends'))
        parents = between.slice('extends'.length).split(',').map((parent) => parent.trim());
      const close = matching(text, open + 1, '{', '}');
      body = text.slice(open + 1, close);
      head.lastIndex = close + 1;
    } else {
      const equals = text.indexOf('=', head.lastIndex);
      let depth = 0;
      let end = equals + 1;
      for (; end < text.length; ++end) {
        const character = text[end];
        if ('{[(<'.includes(character)) ++depth;
        else if ('}])>'.includes(character)) --depth;
        else if (character === ';' && depth === 0) break;
      }
      body = text.slice(equals + 1, end);
      head.lastIndex = end + 1;
    }
    if (found.has(name)) throw new Error(`${name} is declared twice`);
    found.set(name, { kind: keyword === 'type' ? 'alias' : 'interface', extends: parents, body });
  }
  return found;
}

/**
 * Validates runtime payloads against parsed declarations.
 */
export class PayloadChecker {
  /**
   * @param {Map<string, {kind: string, extends: string[], body: string}>} declared -
   *   Declarations from {@link declarations}.
   */
  constructor(declared) {
    this.declared = declared;
    this.parsed = new Map();
    this.reached = new Set();
  }

  /**
   * The type node a declaration resolves to; an interface's members include
   * those it extends.
   * @param {string} name - Declared name.
   * @returns {object} The type node.
   */
  resolve(name) {
    if (this.parsed.has(name)) return this.parsed.get(name);
    const declaration = this.declared.get(name);
    let node;
    if (declaration.kind === 'alias') {
      node = parseType(declaration.body);
    } else {
      const own = parseType(declaration.body, { members: true }).members;
      const inherited = declaration.extends.flatMap((parent) => {
        const base = this.resolve(parent);
        if (base.kind !== 'object') throw new Error(`${name} extends non-object ${parent}`);
        return base.members;
      });
      node = { kind: 'object', members: [...inherited, ...own] };
    }
    this.parsed.set(name, node);
    return node;
  }

  /**
   * Problems with a value against a type node; reached declarations are added
   * to {@link PayloadChecker#reached} only for a value that matches.
   * @param {unknown} value - Runtime value.
   * @param {object} type - Type node.
   * @param {string} path - Location used in problem messages.
   * @returns {string[]} One message per mismatch; empty when the value matches.
   */
  check(value, type, path) {
    const reached = new Set();
    const problems = this.visit(value, type, path, reached);
    if (problems.length === 0) for (const name of reached) this.reached.add(name);
    return problems;
  }

  /**
   * @param {unknown} value - Runtime value.
   * @param {object} type - Type node.
   * @param {string} path - Location used in problem messages.
   * @param {Set<string>} reached - Declarations this match passes through.
   * @returns {string[]} Mismatches.
   */
  visit(value, type, path, reached) {
    const describe = () => (value === null ? 'null' : Array.isArray(value) ? 'array'
      : value?.constructor?.name && typeof value === 'object' ? value.constructor.name : typeof value);
    switch (type.kind) {
      case 'primitive':
        if (type.name === 'unknown') return [];
        if (type.name === 'null') return value === null ? [] : [`${path}: expected null, got ${describe()}`];
        if (type.name === 'object')
          return value !== null && typeof value === 'object' ? [] : [`${path}: expected object, got ${describe()}`];
        return typeof value === type.name ? [] : [`${path}: expected ${type.name}, got ${describe()}`];
      case 'literal':
        return value === type.value ? [] : [`${path}: expected ${JSON.stringify(type.value)}, got ${JSON.stringify(value)}`];
      case 'array':
        if (!Array.isArray(value)) return [`${path}: expected array, got ${describe()}`];
        return value.flatMap((item, index) => this.visit(item, type.element, `${path}[${index}]`, reached));
      case 'tuple':
        if (!Array.isArray(value) || value.length !== type.elements.length)
          return [`${path}: expected ${type.elements.length}-tuple, got ${describe()}`];
        return value.flatMap((item, index) => this.visit(item, type.elements[index], `${path}[${index}]`, reached));
      case 'object': {
        if (value === null || typeof value !== 'object' || Array.isArray(value))
          return [`${path}: expected object, got ${describe()}`];
        const problems = [];
        const names = new Set(type.members.map((member) => member.name));
        for (const key of Object.keys(value))
          if (!names.has(key)) problems.push(`${path}.${key}: field is not declared`);
        for (const member of type.members) {
          if (!(member.name in value)) {
            if (!member.optional) problems.push(`${path}.${member.name}: required field is missing`);
            continue;
          }
          problems.push(...this.visit(value[member.name], member.type, `${path}.${member.name}`, reached));
        }
        return problems;
      }
      case 'union': {
        let closest = null;
        for (const option of type.options) {
          const branch = new Set();
          const problems = this.visit(value, option, path, branch);
          if (problems.length === 0) {
            for (const name of branch) reached.add(name);
            return [];
          }
          const score = problems.length
            + (problems.some((problem) => problem.startsWith(`${path}: expected`)) ? 1e6 : 0);
          if (!closest || score < closest.score) closest = { score, problems };
        }
        return [`${path}: matches no member of the union`, ...closest.problems];
      }
      case 'ref':
        return this.visitRef(value, type, path, reached, describe);
      default:
        throw new Error(`unknown type node ${type.kind}`);
    }
  }

  /**
   * @param {unknown} value - Runtime value.
   * @param {{name: string, args: object[]}} type - Named reference.
   * @param {string} path - Location used in problem messages.
   * @param {Set<string>} reached - Declarations this match passes through.
   * @param {() => string} describe - Names the value's runtime type.
   * @returns {string[]} Mismatches.
   */
  visitRef(value, type, path, reached, describe) {
    if (this.declared.has(type.name)) {
      const lineage = [type.name];
      while (lineage.length) {
        const name = lineage.pop();
        reached.add(name);
        lineage.push(...this.declared.get(name).extends);
      }
      return this.visit(value, this.resolve(type.name), path, reached);
    }
    if (type.name === 'Array' && type.args.length === 1)
      return this.visit(value, { kind: 'array', element: type.args[0] }, path, reached);
    if (type.name === 'Record' && type.args.length === 2) {
      if (value === null || typeof value !== 'object' || Array.isArray(value))
        return [`${path}: expected record, got ${describe()}`];
      return Object.entries(value).flatMap(([key, item]) =>
        this.visit(item, type.args[1], `${path}.${key}`, reached));
    }
    const builtin = globalThis[type.name];
    if (typeof builtin === 'function' && /Array$/.test(type.name))
      return value instanceof builtin ? [] : [`${path}: expected ${type.name}, got ${describe()}`];
    throw new Error(`${path}: type ${type.name} is not declared`);
  }
}

/**
 * Bound functions in Embind's emitted declarations.
 * @param {string} source - holosphere_bindings.d.ts text.
 * @returns {{owner: string, name: string, params: string[], returns: string}[]}
 *   Each function with its owner: `<Class>` for an instance method,
 *   `<Class>.static` for a class function, `module` for a free function.
 */
export function boundFunctions(source) {
  const functions = [];
  const signature = /^\s*(\w+)\((.*)\): (.+);$/;
  let owner = null;
  let embindModule = false;
  for (const line of source.split(/\r?\n/)) {
    const classHead = line.match(/^export interface (\w+) extends ClassHandle \{$/);
    if (classHead) {
      owner = classHead[1];
      continue;
    }
    if (/^interface EmbindModule \{$/.test(line)) {
      embindModule = true;
      owner = 'module';
      continue;
    }
    const nested = embindModule && line.match(/^ {2}(\w+): \{$/);
    if (nested) {
      owner = `${nested[1]}.static`;
      continue;
    }
    if (/^ {2}\};$/.test(line) && embindModule) {
      owner = 'module';
      continue;
    }
    if (/^\}/.test(line)) {
      owner = null;
      embindModule = false;
      continue;
    }
    const match = owner && line.match(signature);
    if (!match || match[1] === 'new') continue;
    const params = match[2] ? match[2].split(/,\s*/).map((param) => param.replace(/^\w+:\s*/, '')) : [];
    functions.push({ owner, name: match[1], params, returns: match[3] });
  }
  return functions;
}

/**
 * Names of the registered payload aliases: the `type <Name>Value = ...`
 * declarations Embind emits for each registered JavaScript payload.
 * @param {string} source - holosphere_bindings.d.ts text.
 * @returns {Set<string>} Alias names.
 */
export function payloadAliases(source) {
  return new Set([...source.matchAll(/^type (\w+Value) = /gm)].map((match) => match[1]));
}

/**
 * Drives every payload-returning bound function and validates each result.
 * @param {object} module - Instantiated engine module.
 * @param {string} bindingsSource - holosphere_bindings.d.ts text.
 * @param {string} typesSource - engine_types.d.ts text.
 * @returns {string[]} Problems; empty when every payload matches.
 */
export function checkEngine(module, bindingsSource, typesSource) {
  const authored = declarations(typesSource);
  const emitted = declarations(bindingsSource);
  const aliases = payloadAliases(bindingsSource);
  const declared = new Map(authored);
  for (const name of aliases) {
    if (declared.has(name)) throw new Error(`${name} is declared by both files`);
    declared.set(name, emitted.get(name));
  }
  const checker = new PayloadChecker(declared);
  const payloadFunctions = boundFunctions(bindingsSource).filter(({ returns }) =>
    tokenize(returns).some((token) => aliases.has(token)));
  const byKey = new Map(payloadFunctions.map((entry) => [`${entry.owner}.${entry.name}`, entry]));
  const sampled = new Set();
  const problems = [];

  const record = (owner, name, value) => {
    const entry = byKey.get(`${owner}.${name}`);
    if (!entry) return;
    sampled.add(`${owner}.${name}`);
    const found = checker.check(value, parseType(entry.returns), `${owner.replace('.static', '')}.${name}()`);
    if (found.length) problems.push(...found.slice(0, 8));
  };
  const call = (target, owner, name, ...args) => {
    const value = target[name](...args);
    record(owner, name, value);
    return value;
  };
  const sampleSimple = (target, owner) => {
    for (const entry of payloadFunctions) {
      if (entry.owner !== owner || !entry.params.every((param) => param === 'number')) continue;
      call(target, owner, entry.name, ...entry.params.map(() => 0.25));
    }
  };

  sampleSimple(module, 'module');
  sampleSimple(module.HolosphereEngine, 'HolosphereEngine.static');

  const engine = new module.HolosphereEngine();
  try {
    for (const { id } of engine.getEffectIdentities()) {
      engine.setEffect(id);
      engine.drawFrame();
      engine.drawFrame();
      sampleSimple(engine, 'HolosphereEngine');
      const chain = engine.getShaderChainBindings();
      if (!chain) continue;
      try {
        sampleSimple(chain, 'ShaderChainBindings');
        const program = chain.getProgram();
        const installed = call(chain, 'ShaderChainBindings', 'setShaderChain', SAMPLE_CHAIN);
        if (installed.status !== module.ChainStatus.OK)
          problems.push(`setShaderChain rejected the sample chain: ${installed.code}`);
        if (chain.setShaderChainParameters(SAMPLE_PARAMETERS) !== module.ParamSetResult.APPLIED)
          problems.push('setShaderChainParameters rejected the sample writes');
        engine.drawFrame();
        sampleSimple(engine, 'HolosphereEngine');
        sampleSimple(chain, 'ShaderChainBindings');
        if (program) call(chain, 'ShaderChainBindings', 'setShaderChain', program);
      } finally {
        chain.delete();
      }
    }
  } finally {
    engine.delete();
  }

  const meshOps = module.MeshOps;
  sampleSimple(meshOps, 'MeshOps.static');
  for (const { name } of meshOps.getRegistry()) {
    call(meshOps, 'MeshOps.static', 'getRecipe', name);
    const mesh = meshOps.fromSolidName(name);
    if (!mesh) {
      problems.push(`MeshOps.fromSolidName(${JSON.stringify(name)}) returned null for a registry solid`);
      continue;
    }
    try {
      sampleSimple(mesh, 'MeshOps');
    } finally {
      mesh.delete();
    }
  }
  meshOps.clearToolingMemory();

  const palettes = new module.PaletteOps();
  try {
    sampleSimple(palettes, 'PaletteOps');
    for (const { recipe } of palettes.effectPresetsV4()) {
      call(palettes, 'PaletteOps', 'compileAndBakeV4', recipe);
      call(palettes, 'PaletteOps', 'inspectV4', recipe);
    }
  } finally {
    palettes.delete();
  }

  for (const key of byKey.keys())
    if (!sampled.has(key)) problems.push(`${key}: returns a payload but no sample call checks it`);
  for (const name of authored.keys())
    if (!LOADER_TYPES.has(name) && !checker.reached.has(name))
      problems.push(`${name}: declared in engine_types.d.ts but no engine payload reaches it`);
  return problems;
}

if (isMain(import.meta.url)) {
  const [modulePath] = process.argv.slice(2);
  if (!modulePath) {
    console.error('usage: engine_payload_types.mjs <holosphere_wasm.js>');
    await exitAfterStderr(2);
  }
  const absolute = isAbsolute(modulePath) ? modulePath : resolve(modulePath);
  const { default: createModule } = await import(pathToFileURL(absolute).href);
  const module = await createModule();
  const problems = checkEngine(module,
    readFileSync(join(dirname(absolute), 'holosphere_bindings.d.ts'), 'utf8'),
    readFileSync(ENGINE_TYPES, 'utf8'));
  if (problems.length) {
    console.error(`engine_payload_types: ${problems.length} payload/declaration mismatches:`);
    for (const problem of problems.slice(0, 40)) console.error(`  ${problem}`);
    if (problems.length > 40) console.error(`  ... ${problems.length - 40} more`);
    await exitAfterStderr(1);
  }
  console.log('engine_payload_types: every engine payload matches engine_types.d.ts and the emitted aliases.');
}
