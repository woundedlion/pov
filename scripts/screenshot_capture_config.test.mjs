import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  CAPTURE_OFFSETS_MS,
  captureIsBlank,
  captureLitFraction,
  DEFAULT_BLANK_FLOOR,
  DEFAULT_CAPTURE_OFFSET_MS,
  minimumLitPixels,
  captureOffsetMs,
} from './screenshot_capture_config.mjs';

test('configured effects resolve to their fixed capture offsets', () => {
  for (const [effect, offset] of Object.entries(CAPTURE_OFFSETS_MS))
    assert.equal(captureOffsetMs(effect), offset, effect);
});

test('unconfigured effects use the default capture offset', () => {
  assert.equal(captureOffsetMs('HankinSolids'), DEFAULT_CAPTURE_OFFSET_MS);
});

test('a global override takes precedence for every effect', () => {
  for (const effect of [...Object.keys(CAPTURE_OFFSETS_MS), 'HankinSolids'])
    assert.equal(captureOffsetMs(effect, 0), 0);
});

test('the default blank floor rejects one pixel below its minimum lit count', () => {
  const minimum = minimumLitPixels();
  assert.ok(minimum > 0);
  assert.equal(captureIsBlank(minimum - 1), true);
  assert.equal(captureLitFraction(minimum - 1) < DEFAULT_BLANK_FLOOR, true);
});

test('the default blank floor accepts its minimum lit count', () => {
  const minimum = minimumLitPixels();
  assert.equal(captureIsBlank(minimum), false);
  assert.equal(captureLitFraction(minimum) >= DEFAULT_BLANK_FLOOR, true);
});
