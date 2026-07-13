// Freshness gate for the docs/screenshots gallery — runs in CI with no browser
// and no rendering (a full re-capture per build is far too slow). It cannot
// detect a screenshot that has gone visually stale relative to the current
// renderer, but it DOES catch the cheap, common form of rot: roster drift. The
// gallery is installed into daydream and served live, so a PNG missing for a
// newly-registered effect (or an orphan PNG for a removed one) ships a broken or
// misleading gallery. Assert the committed PNG set matches the HS_EFFECT_LIST
// roster exactly — the one source the capture script itself iterates.
import { readdir } from 'node:fs/promises';
import { join } from 'node:path';
import { loadEffectRoster, REPO_ROOT } from './effect_roster.mjs';
import { CAPTURE_OFFSETS_MS } from './screenshot_capture_config.mjs';

const SHOTS_DIR = join(REPO_ROOT, 'docs', 'screenshots');

const roster = await loadEffectRoster();
let files;
try {
  files = await readdir(SHOTS_DIR);
} catch (err) {
  if (err.code !== 'ENOENT') throw err;
  files = []; // no gallery dir at all — every roster effect reads as missing below
}
const pngNames = files
  .filter(f => f.endsWith('.png'))
  .map(f => f.slice(0, -'.png'.length));
// Match case-insensitively so the result is identical on the case-insensitive
// Windows dev FS and on Linux CI; a case-only divergence is reported explicitly
// below rather than masked on one FS and double-counted on the other.
const pngByLower = new Map(pngNames.map(p => [p.toLowerCase(), p]));
const rosterLower = new Set(roster.map(e => e.toLowerCase()));

const missing = []; // registered effect, no PNG at all
const caseMismatch = []; // PNG exists but its name differs from the roster only in case
for (const e of roster) {
  const png = pngByLower.get(e.toLowerCase());
  if (png === undefined) missing.push(e);
  else if (png !== e) caseMismatch.push(`${png}.png vs roster '${e}'`);
}
const orphan = pngNames.filter(p => !rosterLower.has(p.toLowerCase())); // PNG, no effect
const invalidOffsets = Object.entries(CAPTURE_OFFSETS_MS)
  .filter(([effect, ms]) => !roster.includes(effect)
    || !Number.isFinite(ms) || ms < 0)
  .map(([effect, ms]) => `${effect}=${ms}`);

if (missing.length || caseMismatch.length || orphan.length || invalidOffsets.length) {
  if (missing.length)
    console.error(`::error::screenshot gallery is missing PNGs for: ${missing.join(', ')}`);
  if (caseMismatch.length)
    console.error(`::error::screenshot PNG names differ from the roster only in case (case-sensitive on Linux CI): ${caseMismatch.join(', ')}`);
  if (orphan.length)
    console.error(`::error::screenshot gallery has orphan PNGs (no such effect): ${orphan.join(', ')}`);
  if (invalidOffsets.length)
    console.error(`::error::screenshot capture offsets are invalid: ${invalidOffsets.join(', ')}`);
  console.error('Regenerate the gallery with: npm run screenshots');
  console.error('(or capture/remove the specific effects above).');
  process.exitCode = 1;
} else {
  console.log(`Screenshot gallery covers all ${roster.length} registered effects (no orphans).`);
}
