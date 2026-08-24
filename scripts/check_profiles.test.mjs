import { test } from 'node:test';
import assert from 'node:assert/strict';
import { mkdtemp, mkdir, readFile, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import {
  profileDirectories,
  PROFILES_DIR,
  reportsIn,
  validateReport,
} from './check_profiles.mjs';

const validReport = `# Example on-device profile — Teensy 4.0 (2026-08-24, **-O3**)

## Setup

Hardware details.

## Frame cadence

Measured cadence.

## Summary ranking

Measured ranking.

## Harness

Harness details.
`;

test('validateReport accepts the checked-in timing report contract', async () => {
  const errors = [];
  const directories = await profileDirectories(PROFILES_DIR);
  let reportCount = 0;
  for (const directory of directories) {
    const reports = await reportsIn(PROFILES_DIR, directory, errors);
    for (const [key, { date, file }] of reports) {
      const report = await readFile(join(PROFILES_DIR, directory, file), 'utf8');
      validateReport(report, directory, key, date, file, errors);
      ++reportCount;
    }
  }
  assert.ok(directories.includes('O3'));
  assert.ok(directories.includes('retired'));
  assert.ok(directories.includes('shipping'));
  assert.ok(reportCount > 0);
  assert.deepEqual(errors, []);
});

test('validateReport rejects a truncated report', () => {
  const errors = [];
  validateReport('# Example on-device profile', 'O3', 'example', '2026-08-24',
    'profile_example_teensy_2026-08-24.md', errors);
  assert.ok(errors.some(error => error.includes('title date')));
  assert.ok(errors.some(error => error.includes('## Frame cadence')));
});

test('validateReport rejects a plausible two-line report', () => {
  const errors = [];
  const report =
    '# Example on-device profile — Teensy 4.0 (2026-08-24, **-O3**)\n\n';
  validateReport(report, 'O3', 'example', '2026-08-24',
    'profile_example_teensy_2026-08-24.md', errors);
  assert.ok(errors.some(error => error.includes('## Setup')));
  assert.ok(errors.some(error => error.includes('## Harness')));
});

test('validateReport binds the heading date to the filename', () => {
  const errors = [];
  validateReport(validReport.replace('2026-08-24', '2026-08-23'), 'O3',
    'example', '2026-08-24', 'profile_example_teensy_2026-08-24.md', errors);
  assert.ok(errors.some(error => error.includes('title date')));
});

test('reportsIn rejects an empty profile set', async t => {
  const root = await mkdtemp(join(tmpdir(), 'holosphere-profiles-'));
  t.after(() => rm(root, { recursive: true, force: true }));
  await mkdir(join(root, 'O3'));
  await writeFile(join(root, 'O3', 'README.md'), '# Global-O3 profiles\n');
  const errors = [];
  const reports = await reportsIn(root, 'O3', errors);
  assert.equal(reports.size, 0);
  assert.deepEqual(errors, ['O3 has no profile reports']);
});

test('profileDirectories derives timing sets from their local indexes', async t => {
  const root = await mkdtemp(join(tmpdir(), 'holosphere-profiles-'));
  t.after(() => rm(root, { recursive: true, force: true }));
  for (const directory of ['shipping', 'O3', 'retired']) {
    await mkdir(join(root, directory));
    await writeFile(join(root, directory, 'README.md'), `# ${directory}\n`);
  }
  await mkdir(join(root, 'memory'));
  await writeFile(join(root, 'memory', 'arena.md'), '# Arena\n');
  assert.deepEqual(await profileDirectories(root), ['O3', 'retired', 'shipping']);
});
