// Pure validation predicates for the headless WASM smoke test.

import { BAKED_CONSTANT_IDS, engineControlNames, fixedDerivedBinding } from './shader_workbench.mjs';
export { BAKED_CONSTANT_IDS, bakedTopologyFields, LIVE_TOPOLOGY_FIELD, engineControlNames } from './shader_workbench.mjs';

/** Fraction of a sub-ceiling stack capacity treated as the creep budget. */
export const STACK_MAX_FILL = 0.75;

/** @returns {string[]} Sparse HyperLattice metadata and write-seam failures. */
export function hyperLatticePatternProblems(observed, results) {
  const problems = [];
  const { pattern, presetIds, accepted, acceptedValue, rejected, rejectedValue } = observed;
  const labels = ['Cubic', 'Octet Truss', 'Shells'];
  const ids = [0, 1, 6];
  const presets = ['cubic-flight', 'cubic-wide-flight', 'hypercube-flight',
    'octet-flight', 'octet-wide-flight', 'octet-4d-flight',
    'shell-flight', 'shell-close-flight', 'shell-4d-flight'];
  if (JSON.stringify(pattern?.options) !== JSON.stringify(labels))
    problems.push('HyperLattice: Pattern labels differ from the shipping patterns');
  if (JSON.stringify(pattern?.optionValues) !== JSON.stringify(ids))
    problems.push('HyperLattice: Pattern optionValues must preserve IDs 0, 1, 6');
  if (JSON.stringify(presetIds) !== JSON.stringify(presets))
    problems.push('HyperLattice: shipping cubic, octet, and shell presets differ');
  if (accepted !== results.APPLIED || acceptedValue !== 6)
    problems.push('HyperLattice: setting Shells ID 6 did not apply');
  if (rejected !== results.INADMISSIBLE || rejectedValue !== 6)
    problems.push('HyperLattice: unavailable Pattern ID 2 did not leave Shells unchanged');
  return problems;
}

/**
 * @param {object} run
 * @param {number} run.frames Frames rendered per effect.
 * @param {Set<string>} run.darkKeys "Name@WxH" passes that lit no pixel.
 * @returns {string[]} One message per all-black pass.
 */
export function darknessProblems({ frames, darkKeys }) {
  return [...darkKeys].map((key) =>
    `${key}: every pixel is zero after ${frames} frame(s): ` +
    `its draw path or the framebuffer view is dead`);
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
    if (d.optionValues !== undefined) {
      const ids = d.optionValues;
      const validMap = Array.isArray(ids) && Array.isArray(d.options) && ids.length > 0
        && ids.length === d.options.length && new Set(ids).size === ids.length
        && ids.every((id) => Number.isInteger(id) && Math.fround(id) === id
          && id >= d.min && id <= d.max)
        && ids.includes(d.min) && ids.includes(d.max);
      if (!validMap)
        problems.push(`param "${d.name}" has invalid optionValues metadata`);
      else if (!ids.includes(d.value))
        problems.push(`param "${d.name}" value ${d.value} is absent from optionValues`);
    }
  }
  return problems;
}

/**
 * Reports unresolved writable promoted ids and invalid compiled-derived values.
 *
 * Writable parameters apply to the compiled effect by control name. Fixed
 * topology and constants are exempt; derived values validate against their
 * source controls. An unresolved writable id refuses the whole apply, so
 * the preview-versus-compiled comparison writes no value at all. Nothing else
 * pins the two vocabularies together: the digests are computed from the
 * document alone, and the value pin in tests/composed_effect/document_values.h
 * maps ids onto parameter families by chain operator, so a label the alias
 * table does not know still lands.
 *
 * @param {object} run
 * @param {{document: string, effect: string, parameterIds: string[], presets?: object[], descriptor?: object}[]} run.documents
 * @param {Map<string, Set<string>>} run.controls Control names per effect id.
 * @param {Set<string>} run.bakedFields From bakedTopologyFields().
 * @returns {string[]} One message per binding or exemption problem.
 */
export function promotedBindingProblems({ documents, controls, bakedFields }) {
  const problems = [];
  const seenConstants = new Set();
  for (const { document, effect, parameterIds, presets = [], descriptor } of documents) {
    const registered = controls.get(effect);
    if (!registered) {
      problems.push(`${document}: the running roster carries no effect "${effect}"`);
      continue;
    }
    for (const parameterId of parameterIds) {
      const authored = presets.filter((preset) =>
        Object.hasOwn(preset.values ?? {}, parameterId));
      const derived = fixedDerivedBinding(descriptor, parameterId, authored[0]?.values ?? {});
      if (derived) {
        if (authored.length === 0 || authored.some(({ values }) =>
          !fixedDerivedBinding(descriptor, parameterId, values).valid) ||
            (derived.sourceId !== null && !engineControlNames(derived.sourceId)
              .some((name) => registered.has(name)))) {
          problems.push(`${document}: "${parameterId}" must match the compiled affine ` +
            'period derived from its source parameters');
        }
        continue;
      }
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
