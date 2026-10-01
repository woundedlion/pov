import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { loadEffectHeaders, stripComments } from './effect_roster.mjs';
import { basename, dirname, relative, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  compilePatternDocuments,
  loadOperatorCatalog,
} from './pattern_documents.mjs';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');

const promotedHeaders = async () => {
  const paths = await loadEffectHeaders();
  const headers = new Map();
  for (const path of paths) {
    const name = relative(resolve(ROOT, 'effects'), path);
    const source = await readFile(path, 'utf8');
    if (!source.includes('DESCRIPTOR_DIGEST')) continue;
    const id = /EFFECT_ID\s*=\s*"([a-z0-9-]+)"/.exec(source);
    assert.ok(id, `${name} carries a digest but no EFFECT_ID`);
    headers.set(id[1],
      { name, types: [basename(path, '.h')] });
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

test('every promoted effect has a document and product group entry', async () => {
  const headers = await promotedHeaders();
  const documents = await compiledDocuments();
  const roster = stripComments(await readFile(resolve(ROOT, 'targets/effects.h'), 'utf8'));
  const group = /^#define HS_SHADER_PRODUCT_GROUP\(X\)(.*)$/m.exec(roster);
  assert.ok(group, 'missing shader product group');
  const names = [...group[1].matchAll(/X\(\s*(\w+)\s*,/g)].map(match => match[1]);
  const promoted = [...headers.values()].flatMap(header => header.types);
  assert.equal(names.length, new Set(names).size, 'duplicate product group effect');
  assert.deepEqual(names.sort(), promoted.sort(), 'product group differs from digest-carrying effects');
  assert.ok(headers.size > 0, 'no promoted effect header carries a digest');

  for (const [effectId, header] of headers) {
    const entry = documents.get(effectId);
    assert.ok(entry, `effects/${header.name} names no pattern document`);
    const dwells = Object.values(entry.compiled.document.preset_bank.choreography.dwell);
    assert.ok(dwells.length > 0, `patterns/${entry.name} has no dwell entries`);
    assert.ok(entry.compiled.document.preset_bank.presets.length === 1 ||
      entry.compiled.document.preset_bank.edges.length > 0,
      `patterns/${entry.name} has no transition edges`);
    for (const dwell of dwells)
      assert.ok(Number.isInteger(dwell) && dwell > 0,
        `patterns/${entry.name} has invalid choreography dwell`);
  }

  for (const [effectId, entry] of documents) {
    assert.ok(headers.has(effectId),
      `patterns/${entry.name} has no promoted header carrying its digests`);
  }
});
