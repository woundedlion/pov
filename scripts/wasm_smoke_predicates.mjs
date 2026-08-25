// The decisions the WASM smoke test fails CI on, separated from the module it
// drives.
//
// wasm_smoke.mjs cannot run without a built WASM module, so nothing gated the
// judgements it makes with one: the all-black expectation for a short frame
// window, the stack creep budget it compares every high-water mark against, and
// the zip between the two embind parameter streams. Each decides whether a real
// defect reds CI, and each is a pure function of numbers and arrays -- kept
// here, free of the module, so wasm_smoke_predicates.test.mjs can gate them
// with no build.

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
