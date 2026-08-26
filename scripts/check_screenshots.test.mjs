import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  invalidCaptureOffsets,
  partitionGallery,
  SHOTS_DIR,
  checkScreenshots,
} from './check_screenshots.mjs';

test('an exact roster/gallery pairing partitions to nothing', () => {
  assert.deepEqual(partitionGallery(['Voronoi', 'Thrusters'],
    ['Thrusters', 'Voronoi']),
  { missing: [], caseMismatch: [], orphan: [] });
});

test('a registered effect with no PNG is missing', () => {
  const { missing, caseMismatch, orphan } =
    partitionGallery(['Voronoi', 'Thrusters'], ['Voronoi']);
  assert.deepEqual(missing, ['Thrusters']);
  assert.deepEqual(caseMismatch, []);
  assert.deepEqual(orphan, []);
});

test('a PNG naming no registered effect is an orphan', () => {
  const { missing, caseMismatch, orphan } =
    partitionGallery(['Voronoi'], ['Voronoi', 'Retired']);
  assert.deepEqual(missing, []);
  assert.deepEqual(caseMismatch, []);
  assert.deepEqual(orphan, ['Retired']);
});

// The divergence that is invisible on the Windows dev FS and fatal on Linux CI:
// it must be its own class, not silently accepted nor double-counted as a
// missing effect plus an orphan PNG.
test('a case-only divergence is reported as a case mismatch', () => {
  const { missing, caseMismatch, orphan } =
    partitionGallery(['IslamicStars'], ['islamicstars']);
  assert.deepEqual(missing, []);
  assert.deepEqual(orphan, []);
  assert.deepEqual(caseMismatch, ["islamicstars.png vs roster 'IslamicStars'"]);
});

// A case-insensitive dev FS cannot hold both spellings, so a stray lowercase
// copy only ever appears on case-sensitive CI.
test('a stray case variant beside the exact PNG is reported', () => {
  const { missing, caseMismatch, orphan } =
    partitionGallery(['Voronoi'], ['voronoi', 'Voronoi']);
  assert.deepEqual(missing, []);
  assert.deepEqual(orphan, []);
  assert.deepEqual(caseMismatch, ["voronoi.png vs roster 'Voronoi'"]);
});

test('capture offsets naming an unregistered effect are rejected', () => {
  assert.deepEqual(invalidCaptureOffsets(['Voronoi'],
    { Voronoi: 30000, Retired: 30000 }), ['Retired=30000']);
});

test('non-finite and negative capture offsets are rejected', () => {
  assert.deepEqual(
    invalidCaptureOffsets(['A', 'B', 'C', 'D'],
      { A: 0, B: -1, C: Number.NaN, D: Number.POSITIVE_INFINITY }),
    ['B=-1', 'C=NaN', 'D=Infinity']);
});

test('the committed gallery satisfies the gate', async () => {
  const result = await checkScreenshots();
  assert.ok(result.rosterSize > 0);
  assert.deepEqual(result.missing, []);
  assert.deepEqual(result.caseMismatch, []);
  assert.deepEqual(result.orphan, []);
  assert.deepEqual(result.invalidOffsets, []);
  assert.deepEqual(result.invalidImages, []);
});

test('a gallery directory that does not exist reads as fully missing', async () => {
  const result = await checkScreenshots(`${SHOTS_DIR}-does-not-exist`);
  assert.equal(result.missing.length, result.rosterSize);
  assert.deepEqual(result.orphan, []);
  assert.deepEqual(result.invalidImages, []);
});
