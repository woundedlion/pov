// Stored gallery geometry. Captures use this exact thumbnail size so a one-pixel
// browser-layout rounding difference cannot produce an invalid gallery asset.
export const GALLERY_WIDTH = 640;
export const GALLERY_HEIGHT = 539;

export const COVERAGE_SAMPLE_WIDTH = 96;
export const COVERAGE_SAMPLE_HEIGHT = 72;
export const COVERAGE_SAMPLE_PIXELS =
  COVERAGE_SAMPLE_WIDTH * COVERAGE_SAMPLE_HEIGHT;
export const DEFAULT_BLANK_FLOOR = 0.01;

export const DEFAULT_CAPTURE_OFFSET_MS = 30_000;

export const CAPTURE_OFFSETS_MS = Object.freeze({
  RingShower: 10_000,
});

/**
 * Returns the fixed post-hydration capture offset for an effect.
 * @param {string} effect Registered effect name.
 * @param {?number} override Global WAIT_MS override, or null when unset.
 * @returns {number} Capture offset in milliseconds.
 */
export function captureOffsetMs(effect, override = null) {
  return override ?? CAPTURE_OFFSETS_MS[effect] ?? DEFAULT_CAPTURE_OFFSET_MS;
}

/**
 * Returns the lit fraction represented by a coverage sample.
 * @param {number} litPixels Number of lit pixels in the coverage sample.
 * @returns {number} Lit fraction in [0, 1].
 */
export function captureLitFraction(litPixels) {
  return litPixels / COVERAGE_SAMPLE_PIXELS;
}

/**
 * Tests a coverage sample against the configured blank floor.
 * @param {number} litPixels Number of lit pixels in the coverage sample.
 * @param {number} blankFloor Minimum accepted lit fraction.
 * @returns {boolean} True when the capture is too sparse to save.
 */
export function captureIsBlank(
  litPixels,
  blankFloor = DEFAULT_BLANK_FLOOR,
) {
  return captureLitFraction(litPixels) < blankFloor;
}

/**
 * Returns the first integer lit-pixel count accepted by a blank floor.
 * @param {number} blankFloor Minimum accepted lit fraction.
 * @returns {number} Minimum accepted lit-pixel count.
 */
export function minimumLitPixels(blankFloor = DEFAULT_BLANK_FLOOR) {
  return Math.ceil(blankFloor * COVERAGE_SAMPLE_PIXELS);
}
