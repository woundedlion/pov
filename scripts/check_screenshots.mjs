// Freshness gate for the docs/screenshots gallery — runs in CI with no browser
// and no rendering (a full re-capture per build is far too slow). It cannot
// detect a screenshot that has gone visually stale relative to the current
// renderer, but it DOES catch the cheap, common forms of rot: roster drift and
// unreadable assets. The gallery is installed into daydream and served live, so
// a PNG missing for a newly-registered effect (or an orphan PNG for a removed
// one) ships a broken or misleading gallery. Assert the committed PNG set
// matches the HS_EFFECT_LIST roster exactly — the one source the capture script
// itself iterates — and that every file is a valid image of the stored size.
import { readdir, readFile } from 'node:fs/promises';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { loadEffectRoster, REPO_ROOT } from './effect_roster.mjs';
import { inspectPng } from './png_probe.mjs';
import {
  CAPTURE_OFFSETS_MS,
  GALLERY_HEIGHT,
  GALLERY_WIDTH,
} from './screenshot_capture_config.mjs';

export const SHOTS_DIR = join(REPO_ROOT, 'docs', 'screenshots');

/**
 * Partitions gallery PNG basenames against the effect roster.
 *
 * Matching is case-insensitive so the result is identical on the
 * case-insensitive Windows dev FS and on Linux CI; a case-only divergence is
 * reported as its own class rather than masked on one FS and counted as both a
 * missing effect and an orphan PNG on the other. Every PNG sharing a roster
 * entry's lowercased name is judged, so a second spelling beside the exact one
 * is reported rather than shadowed.
 *
 * @param {string[]} roster Registered effect names.
 * @param {string[]} pngNames Gallery PNG basenames, without the extension.
 * @returns {{missing: string[], caseMismatch: string[], orphan: string[]}}
 *   Registered effects with no PNG, PNGs whose name differs from the roster
 *   only in case, and PNGs naming no registered effect.
 */
export function partitionGallery(roster, pngNames) {
  const pngsByLower = new Map();
  for (const p of pngNames) {
    const group = pngsByLower.get(p.toLowerCase());
    if (group) group.push(p);
    else pngsByLower.set(p.toLowerCase(), [p]);
  }
  const rosterLower = new Set(roster.map(e => e.toLowerCase()));
  const missing = [];
  const caseMismatch = [];
  for (const e of roster) {
    const pngs = pngsByLower.get(e.toLowerCase());
    if (pngs === undefined) missing.push(e);
    else for (const png of pngs)
      if (png !== e) caseMismatch.push(`${png}.png vs roster '${e}'`);
  }
  const orphan = pngNames.filter(p => !rosterLower.has(p.toLowerCase()));
  return { missing, caseMismatch, orphan };
}

/**
 * Reports the per-effect capture offsets that name no registered effect or are
 * not a finite, non-negative millisecond count.
 *
 * @param {string[]} roster Registered effect names.
 * @param {Object<string, number>} offsets Effect name -> capture offset in ms.
 * @returns {string[]} `effect=value` for each rejected entry.
 */
export function invalidCaptureOffsets(roster, offsets = CAPTURE_OFFSETS_MS) {
  return Object.entries(offsets)
    .filter(([effect, ms]) => !roster.includes(effect)
      || !Number.isFinite(ms) || ms < 0)
    .map(([effect, ms]) => `${effect}=${ms}`);
}

// A name-only check would pass an empty or bit-rotted file, so every gallery
// PNG is validated as a datastream and pinned to the stored gallery geometry.
async function validateImages(shotsDir, pngNames) {
  const invalid = [];
  for (const name of pngNames) {
    try {
      const { width, height } =
        inspectPng(await readFile(join(shotsDir, `${name}.png`)));
      if (width !== GALLERY_WIDTH || height !== GALLERY_HEIGHT)
        invalid.push(`${name}.png is ${width}x${height}, expected ` +
          `${GALLERY_WIDTH}x${GALLERY_HEIGHT}`);
    } catch (err) {
      invalid.push(`${name}.png: ${err.message}`);
    }
  }
  return invalid;
}

export async function checkScreenshots(shotsDir = SHOTS_DIR) {
  const roster = await loadEffectRoster();
  let files;
  try {
    files = await readdir(shotsDir);
  } catch (err) {
    if (err.code !== 'ENOENT') throw err;
    files = []; // no gallery dir at all — every roster effect reads as missing
  }
  const pngNames = files
    .filter(f => f.endsWith('.png'))
    .map(f => f.slice(0, -'.png'.length));
  return {
    rosterSize: roster.length,
    ...partitionGallery(roster, pngNames),
    invalidOffsets: invalidCaptureOffsets(roster),
    invalidImages: await validateImages(shotsDir, pngNames),
  };
}

async function main() {
  const { rosterSize, missing, caseMismatch, orphan, invalidOffsets, invalidImages } =
    await checkScreenshots();
  if (missing.length || caseMismatch.length || orphan.length || invalidOffsets.length
      || invalidImages.length) {
    if (missing.length)
      console.error(`::error::screenshot gallery is missing PNGs for: ${missing.join(', ')}`);
    if (caseMismatch.length)
      console.error(`::error::screenshot PNG names differ from the roster only in case (case-sensitive on Linux CI): ${caseMismatch.join(', ')}`);
    if (orphan.length)
      console.error(`::error::screenshot gallery has orphan PNGs (no such effect): ${orphan.join(', ')}`);
    if (invalidOffsets.length)
      console.error(`::error::screenshot capture offsets are invalid: ${invalidOffsets.join(', ')}`);
    if (invalidImages.length)
      console.error(`::error::screenshot PNGs are not valid ${GALLERY_WIDTH}x${GALLERY_HEIGHT} images: ${invalidImages.join('; ')}`);
    console.error('Regenerate the gallery with: npm run screenshots');
    console.error('(or capture/remove the specific effects above).');
    process.exitCode = 1;
  } else {
    console.log(`Screenshot gallery covers all ${rosterSize} registered effects ` +
      `(no orphans) and every PNG decodes at ${GALLERY_WIDTH}x${GALLERY_HEIGHT}.`);
  }
}

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href)
  await main();
