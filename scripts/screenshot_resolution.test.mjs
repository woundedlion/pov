import { test } from 'node:test';
import assert from 'node:assert/strict';
import { descendToHonoredResolution } from './screenshot_resolution.mjs';

// Returns a load() stub that reports `offers[resolution]` as the selected
// effect, recording the order it was called in.
function stubLoader(offers, calls = []) {
  return async (effect, resolution) => {
    calls.push(resolution);
    return Object.prototype.hasOwnProperty.call(offers, resolution)
      ? offers[resolution] : 'FallbackEffect';
  };
}

test('stops at the first resolution that offers the effect', async () => {
  const calls = [];
  const result = await descendToHonoredResolution(
    'RingSpin', ['big', 'mid', 'small'],
    stubLoader({ big: 'RingSpin', mid: 'RingSpin' }, calls));
  assert.deepEqual(result, { resolution: 'big', honored: true });
  assert.deepEqual(calls, ['big']);
});

test('descends past resolutions the app falls back on', async () => {
  const calls = [];
  const result = await descendToHonoredResolution(
    'RingShower', ['big', 'mid', 'small'],
    stubLoader({ small: 'RingShower' }, calls));
  assert.deepEqual(result, { resolution: 'small', honored: true });
  assert.deepEqual(calls, ['big', 'mid', 'small']);
});

test('reports no honored resolution when every load falls back', async () => {
  const calls = [];
  const result = await descendToHonoredResolution(
    'Ghost', ['big', 'small'], stubLoader({}, calls));
  assert.deepEqual(result, { resolution: 'small', honored: false });
  assert.deepEqual(calls, ['big', 'small']);
});

test('an empty resolution list honors nothing and loads nothing', async () => {
  const calls = [];
  const result = await descendToHonoredResolution(
    'RingSpin', [], stubLoader({}, calls));
  assert.deepEqual(result, { resolution: null, honored: false });
  assert.deepEqual(calls, []);
});

test('a page reporting no effect param is never honored', async () => {
  const result = await descendToHonoredResolution(
    'RingSpin', ['big', 'small'], stubLoader({ big: null, small: undefined }));
  assert.equal(result.honored, false);
});

test('a name the requested effect merely prefixes is not honored', async () => {
  const result = await descendToHonoredResolution(
    'Ring', ['big'], stubLoader({ big: 'RingSpin' }));
  assert.equal(result.honored, false);
});
