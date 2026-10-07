// Preloaded through NODE_OPTIONS: counts the node:assert calls each test
// process runs and writes the count to $HS_ASSERTION_COUNTS.
import assert from 'node:assert';
import { randomUUID } from 'node:crypto';
import { writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { afterEach, beforeEach } from 'node:test';
import { registerHooks, syncBuiltinESMExports } from 'node:module';
import { fileURLToPath } from 'node:url';

const dir = process.env.HS_ASSERTION_COUNTS;
const file = process.argv[1];
if (dir && file && process.env.NODE_TEST_CONTEXT) {
  const loaded = new Set();
  registerHooks({
    load(url, context, nextLoad) {
      if (url.startsWith('file:')) loaded.add(fileURLToPath(url));
      return nextLoad(url, context);
    },
  });
  let count = 0;
  const caseBaselines = [];
  let emptyCases = 0;
  // Wraps every lowercase function property of both specifiers' objects;
  // capitalized keys are classes. A wrapper calls the captured function, so an
  // alias held by two properties scores once. A bare `assert(x)` is bound
  // inside the builtin and scores nothing.
  for (const target of [assert, assert.strict]) {
    for (const key of Object.keys(target)) {
      const fn = target[key];
      if (typeof fn !== 'function' || key === 'strict' || !/^[a-z]/.test(key))
        continue;
      target[key] = (...args) => {
        count += 1;
        return fn.apply(target, args);
      };
    }
  }
  syncBuiltinESMExports();
  beforeEach(() => {
    caseBaselines.push(count);
  });
  afterEach(() => {
    if (count === caseBaselines.pop()) emptyCases += 1;
  });
  // Not the pid: pids recycle within one run and would overwrite a count.
  process.on('exit', () => {
    writeFileSync(
      join(dir, `${randomUUID()}.json`),
      JSON.stringify({ file, count, emptyCases, loaded: [...loaded] }),
    );
  });
}
