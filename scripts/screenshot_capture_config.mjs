export const DEFAULT_CAPTURE_OFFSET_MS = 30_000;

export const CAPTURE_OFFSETS_MS = Object.freeze({
  DistortedRing: 5_000,
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
