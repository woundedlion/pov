import { readdir, readFile } from 'node:fs/promises';
import { join } from 'node:path';
import {
  loadEffectRoster,
  loadPhantasmEffectRoster,
  REPO_ROOT,
} from './effect_roster.mjs';

const PROFILES_DIR = join(REPO_ROOT, 'docs', 'profiles');
const REPORT_RE = /^profile_([a-z0-9]+)_teensy_\d{4}-\d{2}-\d{2}\.md$/;

const errors = [];

async function reportsIn(directory) {
  const entries = await readdir(join(PROFILES_DIR, directory));
  const reports = new Map();
  for (const file of entries) {
    if (file === 'README.md') continue;
    const match = REPORT_RE.exec(file);
    if (!match) {
      errors.push(`${directory}/${file} is not a dated profile report`);
      continue;
    }
    const key = match[1];
    if (reports.has(key))
      errors.push(`${directory} has multiple reports for ${key}`);
    reports.set(key, file);
  }
  return reports;
}

function compareSets(subject, actual, expected) {
  const missing = [...expected].filter(name => !actual.has(name));
  const orphan = [...actual].filter(name => !expected.has(name));
  if (missing.length) errors.push(`${subject} is missing: ${missing.join(', ')}`);
  if (orphan.length) errors.push(`${subject} has orphans: ${orphan.join(', ')}`);
}

function linkedReports(text, prefix) {
  const escaped = prefix.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  const pattern = new RegExp(
    `\\]\\(${escaped}(profile_[a-z0-9]+_teensy_\\d{4}-\\d{2}-\\d{2}\\.md)\\)`,
    'g');
  return new Set([...text.matchAll(pattern)].map(match => match[1]));
}

function checkCount(text, pattern, expected, subject) {
  const match = pattern.exec(text);
  if (!match) errors.push(`${subject} count is not stated`);
  else if (Number(match[1]) !== expected)
    errors.push(`${subject} states ${match[1]}, roster has ${expected}`);
}

const [effectRoster, phantasmRoster] = await Promise.all([
  loadEffectRoster(),
  loadPhantasmEffectRoster(),
]);
const effectKeys = new Set(effectRoster.map(name => name.toLowerCase()));
const phantasmKeys = new Set(phantasmRoster.map(name => name.toLowerCase()));
const shipping = await reportsIn('shipping');
const o3 = await reportsIn('O3');
const retired = await reportsIn('retired');

compareSets('shipping profiles', new Set(shipping.keys()), phantasmKeys);
for (const name of o3.keys()) {
  if (!phantasmKeys.has(name))
    errors.push(`O3 profile names a non-Phantasm effect: ${name}`);
}
for (const name of retired.keys()) {
  if (effectKeys.has(name))
    errors.push(`retired profile still names a registered effect: ${name}`);
}

for (const [directory, reports] of [['shipping', shipping], ['O3', o3],
  ['retired', retired]]) {
  for (const [key, file] of reports) {
    const report = await readFile(join(PROFILES_DIR, directory, file), 'utf8');
    const title = /^# (\w+) on-device profile\b/.exec(report)?.[1];
    if (!title || title.toLowerCase() !== key)
      errors.push(`${directory}/${file} title does not match its filename`);
  }
}

const [index, shippingIndex, o3Index, retiredIndex] = await Promise.all([
  readFile(join(PROFILES_DIR, 'README.md'), 'utf8'),
  readFile(join(PROFILES_DIR, 'shipping', 'README.md'), 'utf8'),
  readFile(join(PROFILES_DIR, 'O3', 'README.md'), 'utf8'),
  readFile(join(PROFILES_DIR, 'retired', 'README.md'), 'utf8'),
]);
compareSets('main shipping index', linkedReports(index, 'shipping/'),
  new Set(shipping.values()));
compareSets('shipping index', linkedReports(shippingIndex, ''),
  new Set(shipping.values()));
compareSets('main O3 index', linkedReports(index, 'O3/'), new Set(o3.values()));
compareSets('O3 index', linkedReports(o3Index, ''), new Set(o3.values()));
compareSets('main retired index', linkedReports(index, 'retired/'),
  new Set(retired.values()));
compareSets('retired index', linkedReports(retiredIndex, ''),
  new Set(retired.values()));
checkCount(index,
  /On-device timing for the \*\*(\d+) effects in the Phantasm image\*\*/,
  phantasmRoster.length, 'main profile index');
checkCount(shippingIndex, /covering the\s+(\d+)\s+effects/, phantasmRoster.length,
  'shipping profile index');

if (errors.length) {
  for (const error of errors) console.error(`::error::${error}`);
  process.exitCode = 1;
} else {
  console.log(`Profile archive covers all ${phantasmRoster.length} Phantasm effects; `
    + `${o3.size} O3 and ${retired.size} retired reports are indexed.`);
}
