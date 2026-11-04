import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { catalogAssignmentFixture } from './generate_catalog_assignment_fixture.mjs';
import { loadOperatorCatalog } from './pattern_documents.mjs';

test('native catalog assignment fixture reproduces the complete current schema', async () => {
  const catalog = await loadOperatorCatalog();
  const generated = catalogAssignmentFixture(catalog);
  const committed = (await readFile(new URL('../tests/composed_effect/generated_catalog_assignments.h', import.meta.url), 'utf8'))
    .replaceAll('\r\n', '\n');
  assert.equal(generated, committed);
  assert.equal([...generated.matchAll(/inline void catalog_assignments_/g)].length, catalog.operators.length);
  assert.match(generated, /params\.terms\[2\]\.amplitude = /);
  assert.match(generated, /params\.lattice_radius = /);
});

test('assignment fixture refuses unclassified or missing operator families', async () => {
  const catalog = await loadOperatorCatalog();
  assert.throws(() => catalogAssignmentFixture({ ...catalog, operators: catalog.operators.slice(1) }), /catalog differ/);
  assert.throws(() => catalogAssignmentFixture({ ...catalog, operators: [
    ...catalog.operators, { id: 'sample.unknown.v1', params: [] },
  ] }), /catalog differ/);
});
