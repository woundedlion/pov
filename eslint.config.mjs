// JavaScript lint rules for first-party tooling: recommended defect rules plus
// checks for counted test assertions. No stylistic rules.
import js from '@eslint/js';
import globals from 'globals';

export default [
  // eslint reads no .gitignore: skip build trees (emitted emscripten glue) and
  // nested checkouts.
  { ignores: ['build*/**', '.worktrees/**', '.hs-pre-commit.*/**', '.doxygen-awesome/**'] },
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
    files: ['**/*.test.mjs'],
    rules: {
      'no-restricted-syntax': ['error', {
        selector: 'CallExpression[callee.type="Identifier"][callee.name="assert"]',
        message:
          'Call node:assert through a property (assert.ok(x)): a bare assert(x) ' +
          'is not counted by scripts/count_assertions.mjs.',
      }, {
        selector: 'ImportDeclaration[source.value=/^(node:)?assert(\\/strict)?$/] > ImportSpecifier',
        message: 'Use the node:assert default or namespace import so assertions are counted.',
      }, {
        selector: 'MemberExpression[property.name="assert"], MemberExpression[property.value="assert"]',
        message: 'Use node:assert methods; test-context assertions are not counted.',
      }],
    },
  },
];
