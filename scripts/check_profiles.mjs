// Validates archive structure and roster coverage; capture freshness is not checked.
import { readdir, readFile } from 'node:fs/promises';
import { basename, dirname, join, relative, resolve } from 'node:path';
import { isMain } from './exit.mjs';
import {
  loadEffectRoster,
  loadPhantasmEffectRoster,
  REPO_ROOT,
} from './effect_roster.mjs';

export const PROFILES_DIR = join(REPO_ROOT, 'docs', 'profiles');
export const REPORT_RE =
  /^profile_([a-z0-9]+)(?:_([a-z0-9]+(?:_[a-z0-9]+)*))?_teensy_(\d{4}-\d{2}-\d{2})\.md$/;
const REQUIRED_SECTIONS = [
  '## Setup',
  '## Frame cadence',
  '## Summary ranking',
  '## Harness',
];

export async function profileDirectories(profilesDir, errors) {
  const entries = await readdir(profilesDir, { withFileTypes: true });
  const directories = [];
  for (const entry of entries) {
    if (!entry.isDirectory()) {
      if (REPORT_RE.test(entry.name))
        errors.push(`${entry.name} is a profile report outside a timing set`);
      continue;
    }
    const files = await readdir(join(profilesDir, entry.name));
    if (files.includes('README.md')) {
      directories.push(entry.name);
      continue;
    }
    // A report set without an index is unreachable and unchecked; a
    // directory holding no report at all is not a timing set.
    if (files.some(file => REPORT_RE.test(file)))
      errors.push(`${entry.name} has profile reports but no README.md index`);
  }
  return directories.sort();
}

export async function reportsIn(profilesDir, directory, errors) {
  const entries = await readdir(join(profilesDir, directory));
  const reports = [];
  const keys = new Set();
  for (const file of entries) {
    if (file === 'README.md') continue;
    const match = REPORT_RE.exec(file);
    if (!match) {
      errors.push(`${directory}/${file} is not a dated profile report`);
      continue;
    }
    const key = match[1];
    const variant = match[2] ?? '';
    const identity = variant ? `${key}_${variant}` : key;
    if (keys.has(identity))
      errors.push(`${directory} has multiple reports for ${identity}`);
    keys.add(identity);
    reports.push({ key, variant, date: match[3], file });
  }
  if (reports.length === 0)
    errors.push(`${directory} has no profile reports`);
  return reports;
}

export function validateReport(report, directory, key, date, file, errors) {
  const subject = `${directory}/${file}`;
  const title =
    /^# ([A-Za-z0-9]+) on-device profile\b[^\r\n]*\((\d{4}-\d{2}-\d{2})(?:,|\))/
      .exec(report);
  if (!title || title[1].toLowerCase() !== key)
    errors.push(`${subject} title does not match its filename`);
  if (!title || title[2] !== date)
    errors.push(`${subject} title date does not match its filename`);

  // Line-anchored: a '### Setup' subheading contains '## Setup' and would
  // otherwise satisfy the requirement and skew every span below.
  const positions = REQUIRED_SECTIONS.map(heading => {
    const found = new RegExp(`^${heading}\\r?$`, 'm').exec(report);
    return found ? found.index : -1;
  });
  for (let i = 0; i < REQUIRED_SECTIONS.length; ++i) {
    if (positions[i] < 0) {
      errors.push(`${subject} is missing ${REQUIRED_SECTIONS[i]}`);
      continue;
    }
    if (i > 0 && positions[i] < positions[i - 1])
      errors.push(`${subject} has out-of-order profile sections`);
    const end = Math.min(report.length,
      ...positions.filter(position => position > positions[i]));
    const bodyStart = report.indexOf('\n', positions[i]) + 1;
    if (bodyStart === 0 || report.slice(bodyStart, end).trim() === '')
      errors.push(`${subject} has no content under ${REQUIRED_SECTIONS[i]}`);
  }
  checkIsrShares(report, subject, errors);
}

const ISR_ROW_RE = /^\s*(?:\|[^|]*?)?(isr_wake|isr_pack|isr_dma_submit)\b/;
const ISR_INLINE_RE =
  /`(isr_wake|isr_pack|isr_dma_submit)`\s+(\d+(?:\.\d+)?)\s*%/g;
const ISR_SUM_PHRASE_RE =
  /ISR counters total|Combined ISR share|ISR CPU shares total|Total (?:measured )?ISR (?:CPU )?share|ISR share totals/i;

/**
 * @brief Rejects prose that adds nested ISR counters into a total.
 * @details `isr_wake` is inclusive of `isr_pack` and `isr_dma_submit`, so a
 *          stated percentage matching wake + pack + submit double-counts them.
 */
export function checkIsrShares(report, subject, errors) {
  const lines = report.replace(/```[^\n]*\n[\s\S]*?```/g, match =>
    match.replace(/[^\n]/g, '')).split(/\r?\n/);
  const groups = [];
  const record = (name, share) => {
    if (name === 'isr_wake') groups.push({ isr_wake: Number(share) });
    else if (groups.length) groups.at(-1)[name] = Number(share);
  };
  for (const line of report.split(/\r?\n/)) {
    const row = ISR_ROW_RE.exec(line);
    if (!row) {
      for (const inline of line.matchAll(ISR_INLINE_RE))
        record(inline[1], inline[2]);
      continue;
    }
    const share = [...line.matchAll(/(\d+(?:\.\d+)?)\s*%/g)].at(-1);
    if (share) record(row[1], share[1]);
  }
  const sums = groups
    .filter(group => 'isr_pack' in group && 'isr_dma_submit' in group)
    .map(group => ({
      wake: group.isr_wake,
      sum: group.isr_wake + group.isr_pack + group.isr_dma_submit,
    }));
  const paragraphs = lines.filter(line => !ISR_ROW_RE.test(line)).join('\n')
    .split(/\n\s*\n/);
  for (const paragraph of paragraphs) {
    if (!/ISR|interrupt|isr_|wake/i.test(paragraph)) continue;
    if (ISR_SUM_PHRASE_RE.test(paragraph)) {
      errors.push(`${subject} sums nested ISR counters; report the inclusive isr_wake share`);
      continue;
    }
    for (const match of paragraph.matchAll(/(\d+\.\d+)(?:\s*[–-]\s*(\d+\.\d+))?\s*%/g)) {
      for (const text of [match[1], match[2]].filter(Boolean)) {
        const value = Number(text);
        const tolerance = 0.5 * 10 ** -text.split('.')[1].length + 0.011;
        if (sums.some(({ wake, sum }) => Math.abs(value - sum) <= tolerance
            && Math.abs(value - wake) > tolerance))
          errors.push(`${subject} states ${text}% = isr_wake + isr_pack + isr_dma_submit; `
            + 'pack and submit are nested in wake, report the inclusive isr_wake share');
      }
    }
  }
}

function compareSets(subject, actual, expected, errors) {
  const missing = [...expected].filter(name => !actual.has(name));
  const orphan = [...actual].filter(name => !expected.has(name));
  if (missing.length) errors.push(`${subject} is missing: ${missing.join(', ')}`);
  if (orphan.length) errors.push(`${subject} has orphans: ${orphan.join(', ')}`);
}

function linkedReports(text, prefix) {
  const escaped = prefix.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  const pattern = new RegExp(
    `\\]\\(${escaped}(${REPORT_RE.source.slice(1, -1)})\\)`,
    'g');
  return new Set([...text.matchAll(pattern)].map(match => match[1]));
}


export async function checkProfileLinks(profilesDir, errors) {
  async function visit(directory) {
    for (const entry of await readdir(directory, { withFileTypes: true })) {
      const path = join(directory, entry.name);
      if (entry.isDirectory()) await visit(path);
      else if (entry.name.endsWith('.md')) {
        const text = (await readFile(path, 'utf8')).replace(/```[^\n]*\n[\s\S]*?```/g, '');
        for (const match of text.matchAll(/\]\(([^\s)]+)(?:\s+"[^"]*")?\)/g)) {
          const link = match[1];
          if (/^[a-z][a-z0-9+.-]*:|^\//i.test(link)) continue;
          const [file, fragment] = link.split('#');
          const target = resolve(dirname(path), decodeURIComponent(file || basename(path)));
          if (relative(profilesDir, target).startsWith('..')) continue;
          let body;
          try { body = await readFile(target, 'utf8'); }
          catch { errors.push(`${relative(profilesDir, path)} has missing link target ${link}`); continue; }
          if (!fragment) continue;
          const anchors = new Set();
          const seen = new Map();
          for (const heading of body.matchAll(/^#{1,6}\s+(.+?)\s*#*$/gm)) {
            const slug = heading[1].toLowerCase().replace(/\[([^\]]+)\]\([^)]*\)/g, '$1')
              .replace(/<[^>]*>/g, '').replace(/[^\p{L}\p{N}_\-\s]/gu, '').replace(/\s/g, '-');
            const count = seen.get(slug) ?? 0;
            seen.set(slug, count + 1);
            anchors.add(count ? `${slug}-${count}` : slug);
          }
          for (const id of body.matchAll(/\b(?:id|name)=["']([^"']+)["']/g)) anchors.add(id[1]);
          if (!anchors.has(decodeURIComponent(fragment)))
            errors.push(`${relative(profilesDir, path)} has missing link anchor ${link}`);
        }
      }
    }
  }
  await visit(profilesDir);
}

function metric(cell) {
  const color = /[🟢🟡🔴]/u.exec(cell)?.[0];
  const values = cell.match(/\d+(?:\.\d+)?/g)?.map(Number);
  return JSON.stringify([color, values]);
}

export function checkIndexCells(index, link, report, subject, errors) {
  const cells = /README cells:\s*peak\s+(.+?),\s*spilled\s+(.+?)\.?(?:\r?\n|$)/.exec(report);
  let headings = [];
  for (const line of index.split(/\r?\n/)) {
    if (!line.startsWith('|')) { headings = []; continue; }
    const row = line.split('|').slice(1, -1).map(cell => cell.trim());
    if (row.some(cell => /peak.*ms/i.test(cell))) { headings = row; continue; }
    if (!line.includes(`](${link})`)) continue;
    const prefix = link.startsWith('O3/') ? 'O3' : link.startsWith('shipping/') ? 'Ship' : '';
    const peak = headings.findIndex(cell => new RegExp(`^${prefix ? prefix + ' ' : ''}peak.*ms$`, 'i').test(cell));
    const spill = headings.findIndex(cell => new RegExp(`^${prefix ? prefix + ' ' : ''}spilled(?:/.*)?$`, 'i').test(cell));
    if (!headings.length) continue;
    if (peak < 0 || spill < 0) {
      errors.push(`${subject} table links ${link} but has no ${prefix ? prefix + ' ' : ''}peak/spilled measurement columns`);
      continue;
    }
    if (!cells) {
      errors.push(`${subject} links ${link} under peak/spill columns but the report has no "README cells:" line`);
      continue;
    }
    if (metric(row[peak]) !== metric(cells[1])) errors.push(`${subject} peak differs from ${link} README cells`);
    if (metric(row[spill]) !== metric(cells[2])) errors.push(`${subject} spill differs from ${link} README cells`);
  }
}

function checkCount(text, pattern, expected, subject, errors) {
  const match = pattern.exec(text);
  if (!match) errors.push(`${subject} count is not stated`);
  else if (Number(match[1]) !== expected)
    errors.push(`${subject} states ${match[1]}, roster has ${expected}`);
}

export async function checkProfiles(profilesDir = PROFILES_DIR) {
  const errors = [];
  const [effectRoster, phantasmRoster, directories, index] = await Promise.all([
    loadEffectRoster(),
    loadPhantasmEffectRoster(),
    profileDirectories(profilesDir, errors),
    readFile(join(profilesDir, 'README.md'), 'utf8'),
  ]);
  await checkProfileLinks(profilesDir, errors);
  const effectKeys = new Set(effectRoster.map(name => name.toLowerCase()));
  const phantasmKeys = new Set(phantasmRoster.map(name => name.toLowerCase()));
  const reportSets = new Map();

  for (const directory of directories)
    reportSets.set(directory, await reportsIn(profilesDir, directory, errors));

  const shipping = reportSets.get('shipping');
  const o3 = reportSets.get('O3');
  const retired = reportSets.get('retired');
  if (!shipping) errors.push('shipping profile directory is missing');
  else {
    compareSets('shipping profiles',
      new Set(shipping.filter(report => !report.variant).map(report => report.key)),
      phantasmKeys, errors);
    for (const { key, variant } of shipping) {
      if (variant && !phantasmKeys.has(key))
        errors.push(`shipping profile names a non-Phantasm effect: ${key}`);
    }
  }
  if (!o3) errors.push('O3 profile directory is missing');
  else {
    for (const { key } of o3) {
      if (!phantasmKeys.has(key))
        errors.push(`O3 profile names a non-Phantasm effect: ${key}`);
    }
  }
  if (retired) {
    for (const { key } of retired) {
      if (effectKeys.has(key))
        errors.push(`retired profile still names a registered effect: ${key}`);
    }
  }

  for (const [directory, reports] of reportSets) {
    for (const { key, date, file } of reports) {
      const report = await readFile(join(profilesDir, directory, file), 'utf8');
      validateReport(report, directory, key, date, file, errors);
      checkIndexCells(index, `${directory}/${file}`, report, 'main profile index', errors);
      const localIndex = await readFile(join(profilesDir, directory, 'README.md'), 'utf8');
      checkIndexCells(localIndex, file, report, `${directory} index`, errors);
    }
    const directoryIndex = await readFile(
      join(profilesDir, directory, 'README.md'), 'utf8');
    const filenames = new Set(reports.map(report => report.file));
    compareSets(`main ${directory} index`,
      linkedReports(index, `${directory}/`), filenames, errors);
    compareSets(`${directory} index`, linkedReports(directoryIndex, ''),
      filenames, errors);
  }

  checkCount(index,
    /On-device timing for the \*\*(\d+) effects in the Phantasm image\*\*/,
    phantasmRoster.length, 'main profile index', errors);
  if (shipping) {
    const shippingIndex = await readFile(
      join(profilesDir, 'shipping', 'README.md'), 'utf8');
    checkCount(shippingIndex, /covering the\s+(\d+)\s+effects/,
      phantasmRoster.length, 'shipping profile index', errors);
  }

  return {
    directories,
    errors,
    o3Count: o3?.length ?? 0,
    phantasmCount: phantasmRoster.length,
    retiredCount: retired?.length ?? 0,
  };
}

async function main() {
  const result = await checkProfiles();
  if (result.errors.length) {
    for (const error of result.errors) console.error(`::error::${error}`);
    process.exitCode = 1;
  } else {
    console.log(
      `Profile archive structure covers all ${result.phantasmCount} Phantasm effects; `
      + `${result.o3Count} O3 and ${result.retiredCount} retired reports are indexed. Capture freshness is not checked.`,
    );
  }
}

if (isMain(import.meta.url))
  await main();
