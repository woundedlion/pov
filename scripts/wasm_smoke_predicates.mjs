// The decisions the WASM smoke test fails CI on, separated from the module it
// drives.
//
// wasm_smoke.mjs cannot run without a built WASM module, so nothing gated the
// judgements it makes with one: the all-black expectation for a short frame
// window, the stack creep budget it compares every high-water mark against, the
// zip between the two embind parameter streams, and the resolution of every
// promoted document parameter id against its effect's controls. Each decides
// whether a real defect reds CI, and each is a pure function of numbers,
// strings and arrays -- kept here, free of the module, so
// wasm_smoke_predicates.test.mjs can gate them with no build.

/** Effect whose rings expand from zero radius, so a short window is black. */
export const DARK_EXEMPT = 'RingShower';

/**
 * Frame counts bracketing the exemption. Below the first, the exempt effect is
 * certainly dark; at or above the second, certainly lit. The per-load RNG
 * reseed moves the crossing by a few frames, so in between either is accepted.
 */
export const DARK_BAND = [24, 48];

/** Fraction of a sub-ceiling stack capacity treated as the creep budget. */
export const STACK_MAX_FILL = 0.75;

/**
 * @param {number} frames Frames rendered per effect.
 * @returns {'exempt-dark'|'either'|'lit'} What the exempt effect must be.
 */
export function darknessExpectation(frames) {
  if (frames < DARK_BAND[0]) return 'exempt-dark';
  return frames >= DARK_BAND[1] ? 'lit' : 'either';
}

/**
 * Every all-black verdict that does not match the expectation.
 *
 * @param {object} run
 * @param {number} run.frames Frames rendered per effect.
 * @param {[number, number][]} run.resolutions Every swept resolution.
 * @param {Set<string>} run.sweptEffects Effect names actually rendered.
 * @param {Set<string>} run.darkKeys "Name@WxH" passes that lit no pixel.
 * @returns {string[]} One message per problem; empty means the sweep agrees.
 */
export function darknessProblems({ frames, resolutions, sweptEffects, darkKeys }) {
  const problems = [];
  // The exemption is only meaningful while it names a live roster entry; a
  // rename would otherwise silently drop the pin.
  if (!sweptEffects.has(DARK_EXEMPT)) {
    problems.push(`the all-black exemption names "${DARK_EXEMPT}", which is not in the ` +
      `rendered roster — the exemption is stale`);
  }
  const expectation = darknessExpectation(frames);
  const wantDark = new Set(expectation === 'exempt-dark'
    ? resolutions.map(([w, h]) => `${DARK_EXEMPT}@${w}x${h}`) : []);
  for (const key of darkKeys) {
    if (wantDark.has(key)) continue;
    if (expectation === 'either' && key.startsWith(`${DARK_EXEMPT}@`)) continue;
    problems.push(`${key}: every pixel is zero after ${frames} frame(s) — ` +
      `its draw path or the framebuffer view is dead`);
  }
  for (const key of wantDark) {
    if (!darkKeys.has(key)) {
      problems.push(`${key}: lit after ${frames} frame(s), but the all-black ` +
        `exemption still claims it is dark`);
    }
  }
  return problems;
}

/**
 * The byte budget a stack high-water mark must stay under.
 *
 * The stack traps nowhere and its mark saturates at capacity, so `hwm >
 * capacity` can never fire; this is the creep tripwire instead. The absolute
 * ceiling is meaningful against any build's stack size, and the capacity
 * fraction covers a hypothetical stack smaller than the ceiling — a degenerate
 * or missing capacity falls back to the ceiling alone, which the caller reports
 * separately.
 *
 * @param {?{capacity: number}} stack The stack region, or a falsy value.
 * @param {number} ceiling Absolute byte ceiling.
 * @param {number} [maxFill] Fraction of capacity allowed.
 * @returns {number}
 */
export function stackCreepBudget(stack, ceiling, maxFill = STACK_MAX_FILL) {
  return Math.min(ceiling,
    stack && stack.capacity > 0 ? stack.capacity * maxFill : ceiling);
}

/**
 * Whether the two embind parameter streams zip.
 *
 * getParameterDefinitions() and getParamValues() share param_marshal.h's
 * ordering and are read in one pass with no drawFrame between, so values[i]
 * must reproduce defs[i].value. Length alone is blind to a transposition, so
 * every index is compared. engine_bindings.h collapses a bool def's value to `raw > 0.5`
 * while the value stream keeps the raw float, so bools are reconstructed rather
 * than compared directly, and they carry no min/max.
 *
 * @param {unknown} defs getParameterDefinitions() result.
 * @param {unknown} values getParamValues() result.
 * @returns {string[]} One message per problem; empty means the seam is intact.
 */
export function paramStreamProblems(defs, values) {
  if (!Array.isArray(defs)) {
    return ['getParameterDefinitions() did not return an array'];
  }
  if (values === null || values === undefined
    || typeof values.length !== 'number') {
    return ['getParamValues() did not return an indexable value stream'];
  }
  const problems = [];
  if (values.length !== defs.length) {
    problems.push(`getParamValues() length ${values.length} != ` +
      `getParameterDefinitions() length ${defs.length} (param order seam drifted)`);
  }
  for (let i = 0; i < defs.length; i++) {
    const d = defs[i];
    if (d === null || typeof d !== 'object') {
      problems.push(`param ${i} is not a definition object`);
      continue;
    }
    if (typeof d.name !== 'string' || d.name.length === 0) {
      problems.push(`param ${i} has no name`);
    }
    const sv = i < values.length ? values[i] : undefined;
    if (i < values.length && !Number.isFinite(sv)) {
      problems.push(`param "${d.name}" value-stream entry ${sv} is not finite`);
    }
    if (typeof d.value === 'boolean') {
      if (i < values.length && d.value !== (sv > 0.5)) {
        problems.push(`param "${d.name}" (index ${i}) def bool ${d.value} ` +
          `!= value-stream ${sv} > 0.5 (param order seam transposed)`);
      }
      continue;
    }
    if (i < values.length && Number.isFinite(sv) && sv !== d.value) {
      problems.push(`param "${d.name}" (index ${i}) def value ${d.value} ` +
        `!= value-stream ${sv} (param order seam transposed)`);
    }
    // Float params carry a finite, ordered range bracketing their value.
    const eps = 1e-4 * (1 + Math.abs(d.max - d.min));
    if (!Number.isFinite(d.min) || !Number.isFinite(d.max) || d.min > d.max) {
      problems.push(`param "${d.name}" has a non-finite/inverted range [${d.min}, ${d.max}]`);
    } else if (!Number.isFinite(d.value) || d.value < d.min - eps || d.value > d.max + eps) {
      problems.push(`param "${d.name}" value ${d.value} outside [${d.min}, ${d.max}]`);
    }
  }
  return problems;
}

/** The one topology field a composed effect leaves live, as a dropdown. */
export const LIVE_TOPOLOGY_FIELD = 'palette-mapping';

/**
 * Field segments a composed effect bakes into its build rather than
 * registering a control for: the catalog's topology parameters, which select
 * an operator's structural variant, less the one field left live.
 *
 * @param {*} catalog The engine operator catalog.
 * @returns {Set<string>} The field segments a fixed apply skips.
 */
export function bakedTopologyFields(catalog) {
  const fields = new Set();
  for (const operator of catalog?.operators ?? []) {
    for (const parameter of operator.params ?? []) {
      if (parameter.topology === true) fields.add(parameter.id);
    }
  }
  fields.delete(LIVE_TOPOLOGY_FIELD);
  return fields;
}

/**
 * Parameter ids a composed effect holds as a compile-time constant: the
 * document carries the value so the chain interpreter reproduces the motion,
 * the compiled build registers no control, and a fixed apply has nothing to
 * write. AshCloud's CAMERA_SPIN_RATE is the only one.
 */
export const BAKED_CONSTANT_IDS = new Set(['camera.spin-speed']);

/** @param {string} value */
const titleWords = (value) => value.split('-')
  .map((part) => (part.length === 0 ? part : part[0].toUpperCase() + part.slice(1)))
  .join(' ');

/**
 * The control names a `<label>.<field>` document parameter id may resolve to
 * on a compiled effect, most specific first.
 *
 * Composed effects register display names off their parameter families, not
 * off document labels, so the two spellings only meet through this table; the
 * simulator's fixed apply path carries the same one.
 *
 * @param {string} parameterId
 * @returns {string[]} The candidate control names.
 */
export function engineControlNames(parameterId) {
  const dot = parameterId.indexOf('.');
  if (dot < 0) return [titleWords(parameterId)];
  const label = parameterId.slice(0, dot);
  const field = parameterId.slice(dot + 1);
  const words = titleWords(field);
  if (label === 'warp1' || label === 'warp2') {
    const slot = `Planar Warp ${label === 'warp1' ? 1 : 2} ${words}`;
    if (['Rotation Rate', 'Translation X', 'Translation Y', 'Scale X', 'Scale Y', 'Shear']
      .includes(words)) return [slot, `Affine ${words}`];
    if (['Radial Scale', 'Radial Phase', 'Angular Phase'].includes(words))
      return [slot, `Polar ${words}`];
    if (['Rotation', 'Cell X', 'Cell Y', 'Offset X', 'Offset Y'].includes(words))
      return [slot, `Mirror ${words}`];
    if (['Strength', 'Frequency', 'Field Angle', 'Scale', 'Vector Angle'].includes(words))
      return [slot, `Warp ${words}`];
    return [slot, words];
  }
  if (label === 'surface') return [`Surface Noise ${words}`];
  if (label === 'camera') return [`Camera ${words}`];
  if (label === 'sample' && field === 'angle-speed') return ['Source Angle Speed'];
  return [words];
}

/**
 * Every promoted document parameter id that names no control on the effect it
 * was promoted to.
 *
 * A promoted document is applied to its compiled effect one parameter at a
 * time, by name; an id that resolves to nothing refuses the whole apply, so
 * the preview-versus-compiled comparison writes no value at all. Nothing else
 * pins the two vocabularies together: the digests are computed from the
 * document alone, and the value pin in tests/test_composed_effect.h maps ids
 * onto parameter families by chain operator, so a label the alias table does
 * not know still lands.
 *
 * @param {object} run
 * @param {{document: string, effect: string, parameterIds: string[]}[]} run.documents
 * @param {Map<string, Set<string>>} run.controls Control names per effect id.
 * @param {Set<string>} run.bakedFields From bakedTopologyFields().
 * @returns {string[]} One message per problem; empty means every id resolves.
 */
export function promotedBindingProblems({ documents, controls, bakedFields }) {
  const problems = [];
  const seenConstants = new Set();
  for (const { document, effect, parameterIds } of documents) {
    const registered = controls.get(effect);
    if (!registered) {
      problems.push(`${document}: the running roster carries no effect "${effect}"`);
      continue;
    }
    for (const parameterId of parameterIds) {
      if (BAKED_CONSTANT_IDS.has(parameterId)) {
        seenConstants.add(parameterId);
        continue;
      }
      if (bakedFields.has(parameterId.slice(parameterId.indexOf('.') + 1))) continue;
      const names = engineControlNames(parameterId);
      if (names.some((name) => registered.has(name))) continue;
      problems.push(`${document}: "${parameterId}" names no control on "${effect}" ` +
        `(tried ${names.map((name) => `"${name}"`).join(', ')}) — applying the ` +
        `document to the compiled build refuses here and writes nothing`);
    }
  }
  // A stale exemption would silently cover an id that has since become
  // registrable, so it only holds while a document still carries it.
  for (const parameterId of BAKED_CONSTANT_IDS) {
    if (!seenConstants.has(parameterId)) {
      problems.push(`the baked-constant exemption names "${parameterId}", which no ` +
        `promoted document carries — the exemption is stale`);
    }
  }
  return problems;
}
