import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, rmSync, writeFileSync } from 'node:fs';
import { dirname, join, relative } from 'node:path';
import { fileURLToPath } from 'node:url';
import { test } from 'node:test';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));
// Under the gitignored build tree: run_tests.mjs resolves the case relative
// to ROOT, and an interrupted run strands nothing `git add -A` would commit.
const FIXTURE_ROOT = join(ROOT, 'build');

test('runner rejects cases without assertions', () => {
  mkdirSync(FIXTURE_ROOT, { recursive: true });
  const fixtureDir = mkdtempSync(join(FIXTURE_ROOT, 'run-tests-fixture-'));
  try {
    const fixture = join(fixtureDir, 'empty.test.mjs');
    writeFileSync(
      fixture,
      `import assert from 'node:assert/strict';
import { test } from 'node:test';

assert.ok(1);
test('empty one', () => {});
test('empty two', () => {});
test('empty three', () => {});
`,
      'utf8',
    );
    const path = relative(ROOT, fixture).replaceAll('\\', '/');
    const run = spawnSync(process.execPath, ['scripts/run_tests.mjs', path], {
      cwd: ROOT,
      encoding: 'utf8',
      env: {
        ...process.env,
        NODE_OPTIONS: '',
        NODE_TEST_CONTEXT: undefined,
        HS_ASSERTION_COUNTS: undefined,
      },
    });

    assert.equal(run.status, 1, run.stdout + run.stderr);
    assert.match(run.stderr, /test cases with no assertions ran in/);
    assert.match(run.stderr, /empty\.test\.mjs/);
  } finally {
    rmSync(fixtureDir, { recursive: true, force: true });
  }
});
