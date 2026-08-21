// Lint config for the repo's Node tooling scripts (scripts/*.mjs), enforced by
// the `lint` job in .github/workflows/ci.yml.
//
// The recommended set only -- the rules that catch defects (undeclared names,
// unreachable code, duplicate keys, unused bindings). No stylistic rules and no
// formatter: the tree passes this unmodified, so the gate reports real breakage
// rather than layout opinions.
import js from '@eslint/js';
import globals from 'globals';

export default [
  // eslint reads no .gitignore, and a wasm build tree holds the emitted
  // emscripten .js glue. Without this, `npm run lint` lints generated code
  // locally while CI (a fresh checkout, no build tree) does not. .worktrees/ is
  // a second full checkout, whose .mjs would otherwise be linted twice.
  { ignores: ['build*/**', '.worktrees/**'] },
  js.configs.recommended,
  {
    files: ['**/*.mjs'],
    languageOptions: {
      ecmaVersion: 'latest',
      sourceType: 'module',
      globals: globals.node,
    },
  },
  {
    // page.evaluate() / addInitScript() bodies are serialized to the browser and
    // run there, so their DOM globals are undeclared in the Node scope.
    files: ['scripts/capture_screenshots.mjs'],
    languageOptions: { globals: { ...globals.node, ...globals.browser } },
  },
  {
    // scripts/count_assertions.mjs cannot wrap node:assert's callable default
    // export, so `assert(x)` is invisible to the nonempty-file check.
    files: ['scripts/*.test.mjs'],
    rules: {
      'no-restricted-syntax': ['error', {
        selector: 'CallExpression[callee.type="Identifier"][callee.name="assert"]',
        message:
          'Call node:assert through a property (assert.ok(x)): a bare assert(x) ' +
          'is not counted by scripts/count_assertions.mjs.',
      }],
    },
  },
];
