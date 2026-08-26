import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readdir, readFile } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  compilePatternDocuments,
  loadOperatorCatalog,
} from './pattern_documents.mjs';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');

const CONSTANT = (name) =>
  new RegExp(`${name}\\s*=\\s*"([0-9a-f]{64})"`);

const promotedHeaders = async () => {
  const names = (await readdir(resolve(ROOT, 'effects')))
    .filter((name) => name.endsWith('.h'));
  const headers = new Map();
  for (const name of names) {
    const source = await readFile(resolve(ROOT, 'effects', name), 'utf8');
    const descriptor = CONSTANT('DESCRIPTOR_DIGEST').exec(source);
    const presetBank = CONSTANT('PRESET_BANK_DIGEST').exec(source);
    if (!descriptor && !presetBank) continue;
    const id = /EFFECT_ID\s*=\s*"([a-z0-9-]+)"/.exec(source);
    assert.ok(id, `${name} carries a digest but no EFFECT_ID`);
    assert.ok(descriptor && presetBank,
      `${name} carries only one of the two digests`);
    headers.set(id[1],
      { name, descriptor: descriptor[1], presetBank: presetBank[1] });
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
  assert.ok(headers.size > 0, 'no promoted effect header carries a digest');

  for (const [effectId, header] of headers) {
    const entry = documents.get(effectId);
    assert.ok(entry, `effects/${header.name} names no pattern document`);
    assert.equal(header.descriptor, entry.compiled.descriptor_digest,
      `effects/${header.name} DESCRIPTOR_DIGEST is stale against patterns/${entry.name}`);
    assert.equal(header.presetBank, entry.compiled.preset_bank_digest,
      `effects/${header.name} PRESET_BANK_DIGEST is stale against patterns/${entry.name}`);
  }

  for (const [effectId, entry] of documents) {
    assert.ok(headers.has(effectId),
      `patterns/${entry.name} has no promoted header carrying its digests`);
  }
});
