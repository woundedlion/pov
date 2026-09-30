// Effect and Phantasm roster readers shared by the documentation gates.
import { readFile, readdir } from 'node:fs/promises';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

export const REPO_ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');

// C++ splices physical lines before removing comments.
export function stripComments(src) {
  return src
    .replace(/\\\r?\n/g, '')
    .replace(/\/\*[\s\S]*?\*\/|\/\/[^\n]*/g,
      (match) => match.startsWith('/*') ? ' ' : '');
}

const WORKBENCH_ONLY_ROWS = new Set([
  'HS_SHADER_WORKBENCH_EFFECT',
  'HS_CHAIN_INTERPRETER_EFFECT',
]);

// Extracts gallery effects from targets/effects.h source text.
export function parseEffectRoster(src) {
  const block = stripComments(src).match(
    /^#define HS_EFFECT_LIST\(X\)(.*)/m);
  if (!block) throw new Error('Could not locate HS_EFFECT_LIST in targets/effects.h');
  const names = [];
  let body = block[1].trim();
  while (body) {
    const row = /^(\w+)\s*\(\s*(\w+)\s*\)/u.exec(body);
    if (!row) throw new Error(`Unknown effect roster row: ${body}`);
    const [, macro, argument] = row;
    if (macro === 'X') names.push(argument);
    else if (argument !== 'X' || !WORKBENCH_ONLY_ROWS.has(macro))
      throw new Error(`Unknown effect roster row: ${row[0]}`);
    body = body.slice(row[0].length).trim();
  }
  if (names.length === 0) throw new Error('HS_EFFECT_LIST parsed to zero effects');
  return names;
}

export async function loadEffectRoster() {
  return parseEffectRoster(
    await readFile(join(REPO_ROOT, 'targets', 'effects.h'), 'utf8'));
}

export function parsePhantasmEffectRoster(src) {
  const block = stripComments(src).match(
    /^#define HS_PHANTASM_EFFECT_LIST\(X\)(.*)/m);
  if (!block)
    throw new Error(
      'Could not locate HS_PHANTASM_EFFECT_LIST in targets/Phantasm/phantasm_playlist.h');
  const names = [...block[1].matchAll(/X\(\s*(\w+)\s*,/g)].map(m => m[1]);
  if (names.length === 0)
    throw new Error('HS_PHANTASM_EFFECT_LIST parsed to zero effects');
  return names;
}

export async function loadPhantasmEffectRoster() {
  return parsePhantasmEffectRoster(await readFile(
    join(REPO_ROOT, 'targets', 'Phantasm', 'phantasm_playlist.h'), 'utf8'));
}

export async function loadEffectHeaders() {
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
  return headers;
}
