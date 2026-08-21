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

// Extracts the X() rows from core/engine/effects.h source text.
export function parseEffectRoster(src) {
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
  // Tolerate whitespace inside the parens: a reformat to `X( Foo )` must not drop
  // rows here, because the same spelling drops them from parseRegisteredEffects too
  // and the cross-check would agree on the truncated roster.
  const names = [...body.matchAll(/X\(\s*(\w+)\s*\)/g)].map(m => m[1]);
  if (names.length === 0) throw new Error('HS_EFFECT_LIST parsed to zero effects');
  return names;
}

export async function loadEffectRoster() {
  return parseEffectRoster(
    await readFile(join(REPO_ROOT, 'core', 'engine', 'effects.h'), 'utf8'));
}

// The standard self-registering effect set, read recursively from effects/.
// The WASM-only workbench surfaces live under workbench/, outside this scan.
// Extracts the REGISTER_EFFECT call sites from one header's source text.
export function parseRegisteredEffects(src) {
  // Strip comments first so a commented-out REGISTER_EFFECT row is not counted
  // (mirrors parseEffectRoster's handling of the X-macro list).
  const stripped = src
    .replace(/\/\*[\s\S]*?\*\//g, '')
    .replace(/\/\/[^\n]*/g, '');
  return [...stripped.matchAll(/REGISTER_EFFECT\(\s*(\w+)\s*\)/g)].map(m => m[1]);
}

export async function loadRegisteredEffects() {
  const dir = join(REPO_ROOT, 'effects');
  const headers = [];
  const visit = async current => {
    for (const entry of await readdir(current, { withFileTypes: true })) {
      const path = join(current, entry.name);
      if (entry.isDirectory()) await visit(path);
      else if (entry.name.endsWith('.h')) headers.push(path);
    }
  };
  await visit(dir);
  const names = new Set();
  for (const path of headers) {
    for (const name of parseRegisteredEffects(await readFile(path, 'utf8'))) {
      names.add(name);
    }
  }
  return [...names];
}
