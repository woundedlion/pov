import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  parseEffectRoster, parsePhantasmEffectRoster,
  loadEffectRoster, loadPhantasmEffectRoster,
  parseEffectIdentities, loadEffectIdentities,
} from './effect_roster.mjs';

test('identity metadata covers gallery classes with canonical IDs and original keys', async () => {
  const identities = await loadEffectIdentities();
  for (const name of await loadEffectRoster())
    assert.ok(identities.some(row => row.name === name));
  assert.deepEqual(identities.find(row => row.name === 'Galaxies'),
    { name: 'Galaxies', id: 'galaxies', seed_key: 'Galaxies' });
});

test('identity parser rejects malformed IDs and cross-effect alias collisions', () => {
  const src = rows => `#define HS_EFFECT_IDENTITY_LIST(X) ${rows}\n`;
  assert.deepEqual(parseEffectIdentities(src('X(Alpha, "alpha", "rng-alpha")')),
    [{ name: 'Alpha', id: 'alpha', seed_key: 'rng-alpha' }]);
  for (const id of ['Alpha', '-alpha', 'alpha-', 'alpha--beta', 'alpha_beta'])
    assert.throws(() => parseEffectIdentities(src(`X(Alpha, "${id}", "rng-alpha")`)), /Invalid canonical/);
  assert.throws(() => parseEffectIdentities(src(
    'X(Alpha, "alpha", "rng-alpha") X(Beta, "alpha", "rng-beta")')), /alias collision/);
  assert.throws(() => parseEffectIdentities(src(
    'X(Alpha, "alpha", "rng-alpha") X(Alpha, "beta", "rng-beta")')), /alias collision/);
  assert.throws(() => parseEffectIdentities(src('')), /zero rows/);
  assert.throws(() => parseEffectIdentities(src('X(Alpha)')), /Unknown effect identity row/);
  assert.throws(() => parseEffectIdentities(''), /Could not locate/);
});

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

test('parsePhantasmEffectRoster returns show rows with and without seconds, in order', () => {
  const src = '/*\n * header\n */\n// show(Commented, 60)\n'
    + 'show(Alpha, 120)\n'
    + '  show( Beta )  // one cycle\n'
    + 'exclude(Gamma)\n';
  assert.deepEqual(parsePhantasmEffectRoster(src), ['Alpha', 'Beta']);
});

test('parsePhantasmEffectRoster throws when there are no show rows', () => {
  assert.throws(() => parsePhantasmEffectRoster('exclude(Alpha)\n'), /has no show rows/);
  assert.throws(() => parsePhantasmEffectRoster(''), /has no show rows/);
});

test('the loaders agree on the checked-in roster', async () => {
  const roster = await loadEffectRoster();
  const phantasm = await loadPhantasmEffectRoster();
  assert.ok(roster.length > 0);
  assert.ok(phantasm.length > 0);
  assert.ok(phantasm.every(name => roster.includes(name)));
});

test('RingShower remains in the full roster and is excluded from Phantasm', async () => {
  assert.ok((await loadEffectRoster()).includes('RingShower'));
  assert.ok(!(await loadPhantasmEffectRoster()).includes('RingShower'));
});


test('roster parsers ignore block-comment openers inside line comments', () => {
  const prefix = '// helpers live in effects/*.h\n';
  assert.deepEqual(parseEffectRoster(prefix + rosterOf('X(Alpha)') + '\n/** doc */'), ['Alpha', 'Last']);
});

test('nested roster macros must identify known workbench-only rows', () => {
  assert.deepEqual(parseEffectRoster(rosterOf('HS_CHAIN_INTERPRETER_EFFECT(X)')), ['Last']);
  assert.throws(() => parseEffectRoster(rosterOf('HS_NEW_EFFECT(X)')), /Unknown effect roster row/);
  assert.throws(() => parseEffectRoster(rosterOf('HS_CHAIN_INTERPRETER_EFFECT(Y)')), /Unknown effect roster row/);
  assert.throws(() => parseEffectRoster(rosterOf('garbage')), /Unknown effect roster row/);
});

test('Phantasm playlist rejects misspelled, unclosed and two-per-line rows', () => {
  for (const row of ['shwo(Alpha, 60)', 'show(Alpha, 60', 'show(Alpha, 60) show(Beta)',
    'show(Alpha, sixty)', 'X(Alpha, 60)'])
    assert.throws(() => parsePhantasmEffectRoster(`show(Ok)\n${row}\n`),
      /Unknown Phantasm playlist row/, row);
});
