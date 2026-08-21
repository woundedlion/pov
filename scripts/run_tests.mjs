// Run `node --test` over the patterns given on the command line and reject any
// test file that executes no cases or no assertions.
import { spawnSync } from 'node:child_process';
import {
  mkdirSync,
  mkdtempSync,
  readFileSync,
  readdirSync,
  rmSync,
} from 'node:fs';
import { tmpdir } from 'node:os';
import { join, relative } from 'node:path';

const COUNTER = new URL('./count_assertions.mjs', import.meta.url).href;
const CASE_REPORTER = new URL('./report_cases.mjs', import.meta.url).href;
const args = process.argv.slice(2);

if (!args.some((arg) => !arg.startsWith('-'))) {
  console.error(
    'run_tests: pass the test file patterns to run, e.g. ' +
      '`node scripts/run_tests.mjs "scripts/*.test.mjs"`.',
  );
  process.exit(1);
}

const suffixes = args
  .filter((arg) => !arg.startsWith('-'))
  .map((pattern) => pattern.slice(pattern.lastIndexOf('*') + 1));
const keyOf = (file) => relative(process.cwd(), file).replaceAll('\\', '/');

const scratch = mkdtempSync(join(tmpdir(), 'holosphere-run-tests-'));
const casesPath = join(scratch, 'cases.json');
const countsDir = join(scratch, 'assertions');
const assertions = new Map();
const cases = new Map();
let status;
let tallied = false;

try {
  mkdirSync(countsDir);
  const run = spawnSync(
    process.execPath,
    [
      '--test',
      '--test-reporter=spec',
      '--test-reporter-destination=stdout',
      `--test-reporter=${CASE_REPORTER}`,
      `--test-reporter-destination=${casesPath}`,
      ...args,
    ],
    {
      stdio: 'inherit',
      env: {
        ...process.env,
        NODE_OPTIONS: `${process.env.NODE_OPTIONS ?? ''} --import ${COUNTER}`,
        HS_ASSERTION_COUNTS: countsDir,
      },
    },
  );
  if (run.error) throw run.error;
  status = run.status ?? 1;

  for (const entry of readdirSync(countsDir)) {
    const { file, count } = JSON.parse(
      readFileSync(join(countsDir, entry), 'utf8'),
    );
    const key = keyOf(file);
    if (key.startsWith('..') || !suffixes.some((suffix) => key.endsWith(suffix)))
      continue;
    assertions.set(key, (assertions.get(key) ?? 0) + count);
  }

  let tallies = '';
  try {
    tallies = readFileSync(casesPath, 'utf8');
  } catch (error) {
    if (error.code !== 'ENOENT') throw error;
  }
  if (tallies.trim() !== '') {
    tallied = true;
    for (const [file, tally] of Object.entries(JSON.parse(tallies)))
      cases.set(keyOf(file), tally.ran);
  }
} finally {
  rmSync(scratch, { recursive: true, force: true });
}

if (status !== 0) process.exit(status);
if (assertions.size === 0) {
  console.error(
    'run_tests: no assertion counts were reported - refusing a green run.',
  );
  process.exit(1);
}
if (!tallied) {
  console.error('run_tests: no case counts were reported - refusing a green run.');
  process.exit(1);
}

const files = new Set([...assertions.keys(), ...cases.keys()]);
const emptyCases = [...files].filter((file) => (cases.get(file) ?? 0) === 0);
const emptyAssertions = [...files].filter(
  (file) => (assertions.get(file) ?? 0) === 0,
);
if (emptyCases.length > 0 || emptyAssertions.length > 0) {
  if (emptyCases.length > 0)
    console.error(
      'run_tests: no test cases ran in:\n' +
        emptyCases
          .sort()
          .map((file) => `  ${file}`)
          .join('\n'),
    );
  if (emptyAssertions.length > 0)
    console.error(
      'run_tests: no assertions ran in:\n' +
        emptyAssertions
          .sort()
          .map((file) => `  ${file}`)
          .join('\n'),
    );
  process.exit(1);
}

const totalCases = [...cases.values()].reduce((sum, count) => sum + count, 0);
const totalAssertions = [...assertions.values()].reduce(
  (sum, count) => sum + count,
  0,
);
console.log(
  `run_tests: ${totalCases} tests passed, ${totalAssertions} assertions ` +
    `across ${files.size} nonempty files.`,
);
