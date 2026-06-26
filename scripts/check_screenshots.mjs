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

const SHOTS_DIR = join(REPO_ROOT, 'docs', 'screenshots');

const roster = await loadEffectRoster();
let files;
try {
  files = await readdir(SHOTS_DIR);
} catch (err) {
  if (err.code !== 'ENOENT') throw err;
  files = []; // no gallery dir at all — every roster effect reads as missing below
}
const pngs = new Set(
  files.filter(f => f.endsWith('.png')).map(f => f.slice(0, -'.png'.length))
);

const missing = roster.filter(e => !pngs.has(e)); // registered effect, no PNG
const orphan = [...pngs].filter(p => !roster.includes(p)); // PNG, no effect

if (missing.length || orphan.length) {
  if (missing.length)
    console.error(`::error::screenshot gallery is missing PNGs for: ${missing.join(', ')}`);
  if (orphan.length)
    console.error(`::error::screenshot gallery has orphan PNGs (no such effect): ${orphan.join(', ')}`);
  console.error('Regenerate the gallery with: npm run screenshots');
  console.error('(or capture/remove the specific effects above).');
  process.exitCode = 1;
} else {
  console.log(`Screenshot gallery covers all ${roster.length} registered effects (no orphans).`);
}
