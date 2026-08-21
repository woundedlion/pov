// A node:test reporter tallying each file's results for scripts/run_tests.mjs.
// TAP attributes no result to a file, so its totals cannot identify an empty
// file.
import { resolve } from 'node:path';

// Reported per file, over `test:pass`/`test:fail`:
//   ran     - tests that reached a result unskipped, which is the file's case
//             count; suites are their tests' parents, not results of their own,
//             so they are not counted.

/**
 * @param {AsyncIterable<{type: string, data: Object}>} source - Runner events.
 * @yields {string} JSON object of `{ ran }` per absolute file path.
 */
export default async function* reportCases(source) {
  const files = new Map();
  for await (const event of source) {
    if (event.type !== 'test:pass' && event.type !== 'test:fail') continue;
    const { file, name, skip, details } = event.data;
    if (!file) continue;
    const tally = files.get(file) ?? { ran: 0 };
    if (
      skip === undefined &&
      details?.type !== 'suite' &&
      resolve(name) !== resolve(file)
    )
      tally.ran += 1;
    files.set(file, tally);
  }
  const tallies = Object.fromEntries(
    [...files].sort(([a], [b]) => (a < b ? -1 : 1)),
  );
  yield `${JSON.stringify(tallies, null, 2)}\n`;
}
