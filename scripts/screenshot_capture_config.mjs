// Stored gallery geometry. The capture downscales to GALLERY_WIDTH and derives
// the height from the live canvas aspect; the CI gate pins the committed result
// to both, so an aspect change lands as a deliberate edit here.
export const GALLERY_WIDTH = 640;
export const GALLERY_HEIGHT = 539;

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
