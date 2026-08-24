import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  CAPTURE_OFFSETS_MS,
  captureIsBlank,
  captureLitFraction,
  DEFAULT_CAPTURE_OFFSET_MS,
  minimumLitPixels,
  captureOffsetMs,
} from './screenshot_capture_config.mjs';

test('configured effects resolve to their fixed capture offsets', () => {
  assert.equal(captureOffsetMs('RingShower'), 10_000);
});

test('unconfigured effects use the default capture offset', () => {
  assert.equal(captureOffsetMs('HankinSolids'), DEFAULT_CAPTURE_OFFSET_MS);
});

test('a global override takes precedence for every effect', () => {
  for (const effect of [...Object.keys(CAPTURE_OFFSETS_MS), 'HankinSolids'])
    assert.equal(captureOffsetMs(effect, 0), 0);
});

test('the default blank floor rejects the last near-empty pixel count', () => {
  assert.equal(minimumLitPixels(), 70);
  assert.equal(captureIsBlank(69), true);
  assert.equal(captureLitFraction(69) < 0.01, true);
});

test('the default blank floor accepts its first boundary pixel count', () => {
  assert.equal(captureIsBlank(70), false);
  assert.equal(captureLitFraction(70) >= 0.01, true);
});
