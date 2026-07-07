// Single source of truth for the effect roster shared by the screenshot tools
// (capture_screenshots.mjs and check_screenshots.mjs): parse the HS_EFFECT_LIST
// X-macro in core/engine/effects.h rather than hand-maintaining a list. That macro is
// the same roster the WASM startup check and the native smoke suite are derived
// from, so the gallery — and its CI freshness gate — can never silently drift
// from the registered effect set when an effect is added or removed.
import { readFile, readdir } from 'node:fs/promises';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

export const REPO_ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');

export async function loadEffectRoster() {
  const src = await readFile(join(REPO_ROOT, 'core', 'effects.h'), 'utf8');
  // Capture the macro body by following its backslash line-continuations rather
  // than relying on a blank line terminating the block (which a reformat could
  // remove). The body runs from `#define HS_EFFECT_LIST(X)` through the last
  // continued line (the first line that does not end in a backslash).
  const block = src.match(/#define HS_EFFECT_LIST\(X\)((?:.*\\\r?\n)*.*)/);
  if (!block) throw new Error('Could not locate HS_EFFECT_LIST in core/engine/effects.h');
  // Strip /* */ block comments then // line comments before extracting names: a
  // commented-out `X(Foo)` row is not in the registered roster, so capturing it
  // would red the gallery run for an otherwise-correct build. Block comments are
  // removed first (and span newlines) so a `/* X(Foo) */` row can't survive into
  // the X() match the way a // row already cannot.
  const body = block[1]
    .replace(/\/\*[\s\S]*?\*\//g, '')
    .replace(/\/\/[^\n]*/g, '');
  const names = [...body.matchAll(/X\((\w+)\)/g)].map(m => m[1]);
  if (names.length === 0) throw new Error('HS_EFFECT_LIST parsed to zero effects');
  return names;
}

// The self-registering effect set, read straight from the headers on disk: every
// effects/*.h carries exactly one REGISTER_EFFECT(ClassName). This is the
// independent witness the HS_EFFECT_LIST roster is cross-checked against — a
// header that registers an effect the X-macro list never names (or vice versa)
// is roster drift the WASM count check cannot see, because an effect never
// #included from core/engine/effects.h is absent from both the registry and the count.
// The REGISTER_EFFECT *macro definition* lives in core/engine/effect_registry.h (not
// effects/), so globbing effects/*.h captures only call sites.
export async function loadRegisteredEffects() {
  const dir = join(REPO_ROOT, 'effects');
  const headers = (await readdir(dir)).filter(f => f.endsWith('.h'));
  const names = new Set();
  for (const f of headers) {
    const src = await readFile(join(dir, f), 'utf8');
    // Strip comments first so a commented-out REGISTER_EFFECT row is not counted
    // (mirrors loadEffectRoster's handling of the X-macro list).
    const stripped = src
      .replace(/\/\*[\s\S]*?\*\//g, '')
      .replace(/\/\/[^\n]*/g, '');
    for (const m of stripped.matchAll(/REGISTER_EFFECT\((\w+)\)/g))
      names.add(m[1]);
  }
  return [...names];
}
