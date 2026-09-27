/* global document */

/** Loads a hydrated effect and confirms selection at its capture offset. */
export async function loadEffectForCapture(page, baseUrl, effect, resolution, offsetMs) {
  const params = new URLSearchParams({ effect, resolution });
  await page.goto(`${baseUrl}?${params.toString()}`,
    { waitUntil: 'load', timeout: 60000 });
  await page.waitForSelector('#canvas', { timeout: 30000 });
  await page.waitForFunction(() => !document.getElementById('loading-overlay'),
    null, { timeout: 60000 });
  const selectedEffect = () => page.evaluate(() =>
    document.querySelector('.effect-button[aria-selected="true"]')?.dataset.effect ?? null);
  const selected = await selectedEffect();
  if (selected !== effect) return selected;
  await page.waitForTimeout(offsetMs);
  return await selectedEffect();
}

/**
 * @typedef {object} DescentResult
 * @property {?string} resolution Last resolution tried, or null if none were.
 * @property {boolean} honored Whether the app selected the requested effect
 *   there. False means the caller must not save a frame.
 */

/**
 * Walks `resolutions` in the given order and stops at the first one whose
 * loaded page reports the requested effect as selected.
 *
 * @param {string} effect Requested effect name.
 * @param {string[]} resolutions Resolutions to try, highest detail first.
 * @param {(effect: string, resolution: string) => Promise<?string>} load
 *   Loads the effect at a resolution and resolves to the effect the app
 *   selected (null when the page reports none).
 * @returns {Promise<DescentResult>}
 */
export async function descendToHonoredResolution(effect, resolutions, load) {
  let resolution = null;
  for (const candidate of resolutions) {
    const selected = await load(effect, candidate);
    resolution = candidate;
    // Strict equality against the request: a null/absent param, a prefix match
    // or the fallback effect all mean this resolution does not offer it.
    if (selected === effect) return { resolution, honored: true };
  }
  return { resolution, honored: false };
}
