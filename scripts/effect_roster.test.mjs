import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  parseEffectRoster, parseRegisteredEffects,
  loadEffectRoster, loadRegisteredEffects,
} from './effect_roster.mjs';

// Both parsers gate CI (check_effect_roster.mjs) and both carry deliberate
// logic — backslash continuation, block comments stripped before line comments,
// and matched whitespace tolerance so a reformat cannot make them agree on a
// truncated roster. Fixtures below cover each of those, plus the failure modes.

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

test('parseEffectRoster drops line-commented rows', () => {
  assert.deepEqual(parseEffectRoster(rosterOf('X(Alpha)', '// X(Dropped)')),
    ['Alpha', 'Last']);
});

test('parseEffectRoster drops block-commented rows, including multi-line ones', () => {
  assert.deepEqual(parseEffectRoster(rosterOf('/* X(Dropped) */ X(Alpha)')),
    ['Alpha', 'Last']);
  assert.deepEqual(parseEffectRoster(rosterOf('/* X(One) \\', 'X(Two) */ X(Alpha)')),
    ['Alpha', 'Last']);
});

// Block comments are stripped first so a `/* ... // ... */` row cannot leave a
// dangling `*/` that resurrects a commented-out X() row.
test('parseEffectRoster strips block comments before line comments', () => {
  assert.deepEqual(parseEffectRoster(rosterOf('/* // X(Dropped)', 'X(AlsoDropped) */ X(Alpha)')),
    ['Alpha', 'Last']);
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

test('parseRegisteredEffects returns every call site in a header', () => {
  assert.deepEqual(
    parseRegisteredEffects('REGISTER_EFFECT(Alpha)\nstruct S {};\nREGISTER_EFFECT(Beta)\n'),
    ['Alpha', 'Beta']);
});

test('parseRegisteredEffects drops commented-out registrations', () => {
  assert.deepEqual(
    parseRegisteredEffects(
      '// REGISTER_EFFECT(Line)\n/* REGISTER_EFFECT(Block) */\nREGISTER_EFFECT(Live)\n'),
    ['Live']);
});

test('parseRegisteredEffects drops multi-line block-commented registrations', () => {
  assert.deepEqual(
    parseRegisteredEffects('/*\nREGISTER_EFFECT(Dropped)\n*/\nREGISTER_EFFECT(Live)\n'),
    ['Live']);
});

// The whitespace tolerance must match parseEffectRoster's: if one parser dropped
// `X( Foo )`-style rows and the other did not, the cross-check would still agree
// on a truncated roster and report nothing.
test('parseRegisteredEffects tolerates the same whitespace parseEffectRoster does', () => {
  assert.deepEqual(parseRegisteredEffects('REGISTER_EFFECT( Alpha )\nREGISTER_EFFECT(\tBeta\t)\n'),
    ['Alpha', 'Beta']);
});

test('parseRegisteredEffects reports none for a header with no registration', () => {
  assert.deepEqual(parseRegisteredEffects('struct NotAnEffect {};\n'), []);
});

// The pure parsers above are only meaningful while they still describe the real
// files the loaders read.
test('the loaders agree on the checked-in roster', async () => {
  const roster = await loadEffectRoster();
  const registered = await loadRegisteredEffects();
  assert.ok(roster.length > 0);
  assert.deepEqual([...roster].sort(), [...registered].sort());
});
