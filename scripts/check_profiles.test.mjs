import { test } from 'node:test';
import assert from 'node:assert/strict';
import { cp, mkdtemp, mkdir, readFile, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import {
  checkProfiles,
  checkProfileLinks,
  checkIndexCells,
  profileDirectories,
  PROFILES_DIR,
  reportsIn,
  validateReport,
} from './check_profiles.mjs';


test('profile links validate supplements, anchors and missing targets', async t => {
  const root = await mkdtemp(join(tmpdir(), 'holosphere-profile-links-'));
  t.after(() => rm(root, { recursive: true, force: true }));
  await mkdir(join(root, 'memory'));
  await writeFile(join(root, 'memory', 'arena.md'), '# Arena\n## Peak usage\n');
  await writeFile(join(root, 'README.md'), '[valid](memory/arena.md#peak-usage)\n');
  const errors = [];
  await checkProfileLinks(root, errors);
  assert.deepEqual(errors, []);
  await writeFile(join(root, 'README.md'), '[bad](memory/arena.md#missing) [gone](hyperlattice_missing.md)\n');
  await checkProfileLinks(root, errors);
  assert.ok(errors.some(error => error.includes('missing link anchor')));
  assert.ok(errors.some(error => error.includes('missing link target')));
});

test('index cells compare paired and local peak and spill values', () => {
  const report = 'README cells: peak 🟢 12.30 (2), spilled 🟢 0/10 (0.00%).\n';
  const index = '| Effect | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled |\n' +
    '| [A](shipping/a.md) / [O3](O3/a.md) | 🟢 12.3 (2) | 🟢 99 | 🟢 0/10 (0%) | 🟢 0/1 (0%) |\n';
  const errors = [];
  checkIndexCells(index, 'shipping/a.md', report, 'main', errors);
  assert.deepEqual(errors, []);
  checkIndexCells(index, 'O3/a.md', report, 'main', errors);
  assert.equal(errors.length, 2);
  const local = '| Effect | Peak ms | Spilled |\n| [A](a.md) | 🟢 12.3 (2) | 🟢 1/10 (10%) |\n';
  checkIndexCells(local, 'a.md', report, 'shipping', errors);
  assert.ok(errors.at(-1).includes('spill differs'));
});

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
  const directories = await profileDirectories(PROFILES_DIR, errors);
  let reportCount = 0;
  for (const directory of directories) {
    const reports = await reportsIn(PROFILES_DIR, directory, errors);
    for (const { key, date, file } of reports) {
      const report = await readFile(join(PROFILES_DIR, directory, file), 'utf8');
      validateReport(report, directory, key, date, file, errors);
      ++reportCount;
    }
  }
  assert.ok(directories.includes('O3'));
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

test('validateReport does not accept a subheading as a required section', () => {
  const errors = [];
  validateReport(validReport.replace('## Setup', '### Setup'), 'O3',
    'example', '2026-08-24', 'profile_example_teensy_2026-08-24.md', errors);
  assert.ok(errors.some(error => error.includes('missing ## Setup')));
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
  assert.equal(reports.length, 0);
  assert.deepEqual(errors, ['O3 has no profile reports']);
});

test('reportsIn keeps both reports when two share an effect key', async t => {
  const root = await mkdtemp(join(tmpdir(), 'holosphere-profiles-'));
  t.after(() => rm(root, { recursive: true, force: true }));
  await mkdir(join(root, 'O3'));
  await writeFile(join(root, 'O3', 'README.md'), '# Global-O3 profiles\n');
  for (const date of ['2026-08-23', '2026-08-24']) {
    await writeFile(join(root, 'O3', `profile_example_teensy_${date}.md`),
      validReport);
  }
  const errors = [];
  const reports = await reportsIn(root, 'O3', errors);
  assert.deepEqual(reports.map(report => report.file), [
    'profile_example_teensy_2026-08-23.md',
    'profile_example_teensy_2026-08-24.md',
  ]);
  assert.deepEqual(errors, ['O3 has multiple reports for example']);
});

test('reportsIn separates supplemental variants and rejects duplicate variants', async t => {
  const root = await mkdtemp(join(tmpdir(), 'holosphere-profiles-'));
  t.after(() => rm(root, { recursive: true, force: true }));
  await mkdir(join(root, 'O3'));
  for (const name of [
    'example_teensy_2026-08-24',
    'example_preset3_teensy_2026-08-24',
    'example_octet_preset5_teensy_2026-08-23',
    'example_octet_preset5_teensy_2026-08-24',
  ]) {
    await writeFile(join(root, 'O3', `profile_${name}.md`), validReport);
  }
  const errors = [];
  const reports = await reportsIn(root, 'O3', errors);
  assert.equal(reports.length, 4);
  assert.ok(reports.every(report => report.key === 'example'));
  assert.deepEqual(new Set(reports.map(report => report.variant)),
    new Set(['', 'preset3', 'octet_preset5']));
  assert.deepEqual(errors, ['O3 has multiple reports for example_octet_preset5']);
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
  const errors = [];
  assert.deepEqual(await profileDirectories(root, errors),
    ['O3', 'retired', 'shipping']);
  assert.deepEqual(errors, []);
});

test('profileDirectories reports a report loose in the archive root', async t => {
  const root = await mkdtemp(join(tmpdir(), 'holosphere-profiles-'));
  t.after(() => rm(root, { recursive: true, force: true }));
  await writeFile(join(root, 'profile_example_teensy_2026-08-24.md'), validReport);
  await writeFile(join(root, 'notes.md'), '# Notes\n');
  const errors = [];
  assert.deepEqual(await profileDirectories(root, errors), []);
  assert.deepEqual(errors,
    ['profile_example_teensy_2026-08-24.md is a profile report outside a timing set']);
});

test('profileDirectories reports an unindexed set of reports', async t => {
  const root = await mkdtemp(join(tmpdir(), 'holosphere-profiles-'));
  t.after(() => rm(root, { recursive: true, force: true }));
  await mkdir(join(root, 'orphan'));
  await writeFile(join(root, 'orphan', 'profile_example_teensy_2026-08-24.md'),
    validReport);
  const errors = [];
  assert.deepEqual(await profileDirectories(root, errors), []);
  assert.deepEqual(errors,
    ['orphan has profile reports but no README.md index']);
});

test('checkProfiles accepts current reports without retired archives', async t => {
  const root = await mkdtemp(join(tmpdir(), 'holosphere-profile-gate-'));
  t.after(() => rm(root, { recursive: true, force: true }));
  await cp(PROFILES_DIR, root, { recursive: true });
  await rm(join(root, 'retired'), { recursive: true, force: true });
  const result = await checkProfiles(root);
  assert.deepEqual(result.errors, []);
  assert.equal(result.retiredCount, 0);
});

test('checkProfiles validates cross-roster and index contracts', async t => {
  assert.deepEqual((await checkProfiles()).errors, []);
  const reports = await reportsIn(PROFILES_DIR, 'shipping', []);
  const first = reports[0];
  const supplement = reports.find(report => report.variant);
  const canonical = reports.find(report =>
    report.key === supplement.key && !report.variant);
  assert.ok(canonical);
  const cases = [
    ['missing shipping directory', async root => {
      await rm(join(root, 'shipping'), { recursive: true });
    }, 'shipping profile directory is missing'],
    ['missing O3 directory', async root => {
      await rm(join(root, 'O3'), { recursive: true });
    }, 'O3 profile directory is missing'],
    ['missing shipping report', async root => {
      await rm(join(root, 'shipping', first.file));
    }, 'shipping profiles is missing: ' + first.key],
    ['supplement cannot replace canonical shipping report', async root => {
      await rm(join(root, 'shipping', canonical.file));
    }, 'shipping profiles is missing: ' + canonical.key],
    ['supplement names an unknown effect', async root => {
      await writeFile(join(root, 'shipping',
        'profile_example_preset3_teensy_2026-08-24.md'), validReport);
    }, 'shipping profile names a non-Phantasm effect: example'],
    ['supplement has a wrong title', async root => {
      await writeFile(join(root, 'shipping', supplement.file), validReport);
    }, 'shipping/' + supplement.file + ' title does not match its filename'],
    ['supplement missing main index link', async root => {
      const path = join(root, 'README.md');
      const text = await readFile(path, 'utf8');
      await writeFile(path, text.replaceAll('shipping/' + supplement.file, 'removed.md'));
    }, 'main shipping index is missing: ' + supplement.file],
    ['supplement missing local index link', async root => {
      const path = join(root, 'shipping', 'README.md');
      const text = await readFile(path, 'utf8');
      await writeFile(path, text.replaceAll(supplement.file, 'removed.md'));
    }, 'shipping index is missing: ' + supplement.file],
    ['orphan shipping report', async root => {
      await writeFile(join(root, 'shipping',
        'profile_example_teensy_2026-08-24.md'), validReport);
    }, 'shipping profiles has orphans: example'],
    ['non-Phantasm reference', async root => {
      await writeFile(join(root, 'O3',
        'profile_example_teensy_2026-08-24.md'), validReport);
    }, 'O3 profile names a non-Phantasm effect: example'],
    ['registered retired effect', async root => {
      await mkdir(join(root, 'retired'), { recursive: true });
      await writeFile(join(root, 'retired', 'README.md'), '# Retired profiles\n');
      await cp(join(root, 'shipping', first.file), join(root, 'retired', first.file));
    }, 'retired profile still names a registered effect: ' + first.key],
    ['missing main index link', async root => {
      const path = join(root, 'README.md');
      const text = await readFile(path, 'utf8');
      await writeFile(path, text.replaceAll('shipping/' + first.file, 'removed.md'));
    }, 'main shipping index is missing: ' + first.file],
    ['missing local index link', async root => {
      const path = join(root, 'shipping', 'README.md');
      const text = await readFile(path, 'utf8');
      await writeFile(path, text.replaceAll(first.file, 'removed.md'));
    }, 'shipping index is missing: ' + first.file],
    ['main count', async root => {
      const path = join(root, 'README.md');
      const text = await readFile(path, 'utf8');
      await writeFile(path, text.replace(/\*\*\d+ effects in the Phantasm image/,
        '**999 effects in the Phantasm image'));
    }, 'main profile index states 999'],
    ['shipping count', async root => {
      const path = join(root, 'shipping', 'README.md');
      const text = await readFile(path, 'utf8');
      await writeFile(path, text.replace(/covering the\s+\d+\s+effects/,
        'covering the 999 effects'));
    }, 'shipping profile index states 999'],
  ];
  for (const [name, mutate, expected] of cases) {
    await t.test(name, async subtest => {
      const root = await mkdtemp(join(tmpdir(), 'holosphere-profile-gate-'));
      subtest.after(() => rm(root, { recursive: true, force: true }));
      await cp(PROFILES_DIR, root, { recursive: true });
      await mutate(root);
      const result = await checkProfiles(root);
      assert.ok(result.errors.some(error => error.startsWith(expected)),
        JSON.stringify(result.errors));
    });
  }
});
