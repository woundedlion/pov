import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  parseEffectRoster, parsePhantasmEffectRoster,
  loadEffectRoster, loadPhantasmEffectRoster,
} from './effect_roster.mjs';

const rosterOf = (...rows) =>
  `#define HS_EFFECT_LIST(X) \\\n${rows.map(r => `  ${r} \\`).join('\n')}\n  X(Last)\n\nint after;\n`;

test('parseEffectRoster follows backslash continuations to the last row', () => {
  assert.deepEqual(parseEffectRoster(rosterOf('X(Alpha)', 'X(Beta)')),
    ['Alpha', 'Beta', 'Last']);
});

test('parseEffectRoster stops at the first uncontinued line', () => {
  const src = '#define HS_EFFECT_LIST(X) \\\n  X(Alpha)\n#define OTHER X(NotAnEffect)\n';
  assert.deepEqual(parseEffectRoster(src), ['Alpha']);
});

test('parseEffectRoster survives CRLF continuations', () => {
  const src = '#define HS_EFFECT_LIST(X) \\\r\n  X(Alpha) \\\r\n  X(Beta)\r\n';
  assert.deepEqual(parseEffectRoster(src), ['Alpha', 'Beta']);
});

test('parseEffectRoster splices continuations before dropping line comments', () => {
  assert.deepEqual(parseEffectRoster(rosterOf('X(Alpha)', '// X(Dropped)')),
    ['Alpha']);
});

test('parseEffectRoster drops block-commented rows, including multi-line ones', () => {
  assert.deepEqual(parseEffectRoster(rosterOf('/* X(Dropped) */ X(Alpha)')),
    ['Alpha', 'Last']);
  assert.deepEqual(parseEffectRoster(rosterOf('/* X(One) \\', 'X(Two) */ X(Alpha)')),
    ['Alpha', 'Last']);
});

test('parseEffectRoster respects line-comment openers inside block comments', () => {
  assert.deepEqual(parseEffectRoster(rosterOf('/* // X(Dropped)', 'X(AlsoDropped) */ X(Alpha)')),
    ['Alpha', 'Last']);
});

// The macro is named in prose all over the tree; a comment spelling out its
// definition must not be read as the definition.
test('parseEffectRoster ignores a commented-out definition of the macro', () => {
  const line = '// #define HS_EFFECT_LIST(X) names every effect\n'
    + rosterOf('X(Alpha)');
  assert.deepEqual(parseEffectRoster(line), ['Alpha', 'Last']);
  const block = '/*\n#define HS_EFFECT_LIST(X) names every effect\n*/\n'
    + rosterOf('X(Alpha)');
  assert.deepEqual(parseEffectRoster(block), ['Alpha', 'Last']);
});

test('parseEffectRoster tolerates whitespace inside the parens', () => {
  assert.deepEqual(parseEffectRoster(rosterOf('X( Alpha )', 'X(\tBeta\t)')),
    ['Alpha', 'Beta', 'Last']);
});

test('parseEffectRoster throws when the macro is missing', () => {
  assert.throws(() => parseEffectRoster('X(Alpha)\n'), /Could not locate HS_EFFECT_LIST/);
});

test('parseEffectRoster throws on an empty roster rather than reporting none', () => {
  assert.throws(() => parseEffectRoster('#define HS_EFFECT_LIST(X)\n'),
    /parsed to zero effects/);
});

test('parsePhantasmEffectRoster accepts literal and derived durations', () => {
  const src = '#define HS_PHANTASM_EFFECT_LIST(X) \\\n'
    + '  X(Alpha, 120) \\\n'
    + '  X(Beta, \\\n'
    + '    hs_preset_window_seconds<Beta>())\n';
  assert.deepEqual(parsePhantasmEffectRoster(src), ['Alpha', 'Beta']);
});

test('parsePhantasmEffectRoster throws when the macro is missing or empty', () => {
  assert.throws(() => parsePhantasmEffectRoster('X(Alpha, 120)\n'),
    /Could not locate HS_PHANTASM_EFFECT_LIST/);
  assert.throws(() => parsePhantasmEffectRoster(
    '#define HS_PHANTASM_EFFECT_LIST(X)\n'), /parsed to zero effects/);
});

// The pure parsers above are only meaningful while they still describe the real
// files the loaders read.
test('the loaders agree on the checked-in roster', async () => {
  const roster = await loadEffectRoster();
  const phantasm = await loadPhantasmEffectRoster();
  assert.ok(roster.length > 0);
  assert.ok(phantasm.length > 0);
  assert.ok(phantasm.every(name => roster.includes(name)));
});

test('roster parsers ignore block-comment openers inside line comments', () => {
  const prefix = '// helpers live in effects/*.h\n';
  assert.deepEqual(parseEffectRoster(prefix + rosterOf('X(Alpha)') + '\n/** doc */'), ['Alpha', 'Last']);
});
