import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { loadEffectHeaders, parseRegisteredEffects } from './effect_roster.mjs';
import { dirname, relative, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  compilePatternDocuments,
  loadOperatorCatalog,
} from './pattern_documents.mjs';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');

const CONSTANT = (name) =>
  new RegExp(`${name}\\s*=\\s*"([0-9a-f]{64})"`);
const INTEGER_CONSTANT = (name) =>
  new RegExp(`${name}\\s*=\\s*(\\d+)`);

const promotedHeaders = async () => {
  const paths = await loadEffectHeaders();
  const headers = new Map();
  for (const path of paths) {
    const name = relative(resolve(ROOT, 'effects'), path);
    const source = await readFile(path, 'utf8');
    const descriptor = CONSTANT('DESCRIPTOR_DIGEST').exec(source);
    const presetBank = CONSTANT('PRESET_BANK_DIGEST').exec(source);
    if (!descriptor && !presetBank) continue;
    const id = /EFFECT_ID\s*=\s*"([a-z0-9-]+)"/.exec(source);
    assert.ok(id, `${name} carries a digest but no EFFECT_ID`);
    assert.ok(descriptor && presetBank,
      `${name} carries only one of the two digests`);
    const dwell = INTEGER_CONSTANT('PRESET_DWELL_FRAMES').exec(source);
    assert.ok(dwell, `${name} carries digests but no PRESET_DWELL_FRAMES`);
    headers.set(id[1],
      { name, types: parseRegisteredEffects(source), descriptor: descriptor[1], presetBank: presetBank[1], dwell: Number(dwell[1]) });
  }
  return headers;
};

const compiledDocuments = async () => {
  const documents = new Map();
  const catalog = await loadOperatorCatalog();
  for (const { name, compiled } of await compilePatternDocuments(catalog)) {
    assert.equal(compiled.status, 'VALID', `patterns/${name} does not compile`);
    // A study document names no effect: only a promoted one pins a header.
    if (typeof compiled.document.effect_id === 'string')
      documents.set(compiled.document.effect_id, { name, compiled });
  }
  return documents;
};

test('every promoted header digest matches its pattern document', async () => {
  const headers = await promotedHeaders();
  const documents = await compiledDocuments();
  const roster = (await readFile(resolve(ROOT, 'targets/effects.h'), 'utf8'))
    .replace(/\\\r?\n/g, '')
    .replace(/\/\*[\s\S]*?\*\//g, ' ')
    .replace(/\/\/[^\n]*/g, '');
  const group = /^#define HS_SHADER_PRODUCT_GROUP\(X\)(.*)$/m.exec(roster);
  assert.ok(group, 'missing shader product group');
  const names = [...group[1].matchAll(/X\(\s*(\w+)\s*,/g)].map(match => match[1]);
  const promoted = [...headers.values()].flatMap(header => header.types);
  assert.equal(names.length, new Set(names).size, 'duplicate product group effect');
  assert.deepEqual(names.sort(), promoted.sort(), 'product group differs from digest-carrying effects');
  const composedEffect = await readFile(
    resolve(ROOT, 'core/render/pullback/composed_effect.h'), 'utf8');
  const segue = /PRESET_SEGUE\s*\{\s*(\d+)/.exec(composedEffect);
  assert.ok(segue, 'composed_effect.h carries no PRESET_SEGUE duration');
  const segueDuration = Number(segue[1]);
  assert.ok(headers.size > 0, 'no promoted effect header carries a digest');

  for (const [effectId, header] of headers) {
    const entry = documents.get(effectId);
    assert.ok(entry, `effects/${header.name} names no pattern document`);
    assert.equal(header.descriptor, entry.compiled.descriptor_digest,
      `effects/${header.name} DESCRIPTOR_DIGEST is stale against patterns/${entry.name}`);
    assert.equal(header.presetBank, entry.compiled.preset_bank_digest,
      `effects/${header.name} PRESET_BANK_DIGEST is stale against patterns/${entry.name}`);
    const dwells = Object.values(entry.compiled.document.preset_bank.choreography.dwell);
    assert.ok(dwells.length > 0, `patterns/${entry.name} has no dwell entries`);
    assert.ok(entry.compiled.document.preset_bank.presets.length === 1 ||
      entry.compiled.document.preset_bank.edges.length > 0,
      `patterns/${entry.name} has no transition edges`);
    for (const dwell of dwells)
      assert.equal(dwell, header.dwell,
        `patterns/${entry.name} choreography dwell differs from effects/${header.name}`);
    for (const edge of entry.compiled.document.preset_bank.edges)
      assert.equal(edge.duration, segueDuration,
        `patterns/${entry.name} edge duration differs from PRESET_SEGUE`);
  }

  for (const [effectId, entry] of documents) {
    assert.ok(headers.has(effectId),
      `patterns/${entry.name} has no promoted header carrying its digests`);
  }
});
