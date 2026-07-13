import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  CAPTURE_OFFSETS_MS,
  DEFAULT_CAPTURE_OFFSET_MS,
  captureOffsetMs,
} from './screenshot_capture_config.mjs';

test('configured effects resolve to their fixed capture offsets', () => {
  assert.equal(captureOffsetMs('DistortedRing'), 5_000);
  assert.equal(captureOffsetMs('RingShower'), 10_000);
});

test('unconfigured effects use the default capture offset', () => {
  assert.equal(captureOffsetMs('HankinSolids'), DEFAULT_CAPTURE_OFFSET_MS);
});

test('a global override takes precedence for every effect', () => {
  for (const effect of [...Object.keys(CAPTURE_OFFSETS_MS), 'HankinSolids'])
    assert.equal(captureOffsetMs(effect, 0), 0);
});
