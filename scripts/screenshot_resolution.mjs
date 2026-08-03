// Resolution descent for the gallery capture: the decision that keeps a
// screenshot from being written under the wrong effect's filename.
//
// The simulator rewrites the URL's effect param to whatever it actually
// selected, so an effect the requested resolution does not offer comes back as
// the app's fallback. Capturing that frame would overwrite a correct PNG with a
// thumbnail of a different effect, which the freshness gate cannot detect (it
// validates names and dimensions, not content). Kept here, free of Playwright,
// so screenshot_resolution.test.mjs can gate it without a browser.

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
