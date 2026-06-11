// Headless Playwright script that loads the WebGL simulator for each effect,
// lets it animate, and saves a PNG screenshot to docs/screenshots/.
// Effects can be overridden via CLI args; otherwise the full EFFECTS list runs.
//
// Usage (from the Holosphere repo root):
//   1. Serve the sibling daydream checkout (see README §"Running the Simulator"):
//          cd ../daydream && python3 -m http.server 8080
//   2. Install the Playwright browser once:  npx playwright install chromium
//   3. Capture the gallery:                  npm run screenshots
//      (equivalently:  node scripts/capture_screenshots.mjs [Effect ...])
//
// SIM_URL overrides the simulator origin (defaults to the README's local
// http.server port); WAIT_MS overrides the per-effect animation settle time.
import { chromium } from 'playwright';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const BASE_URL = process.env.SIM_URL || 'http://localhost:8080/';
const RESOLUTION = 'Phantasm (144x288)';
const OUT_DIR = 'docs/screenshots';
const WAIT_MS = parseInt(process.env.WAIT_MS || '30000', 10);

// Single source of truth for the effect roster: parse the HS_EFFECT_LIST X-macro
// in core/effects.h rather than hand-maintaining a list here. That macro is the
// same roster the WASM startup check and the native smoke suite are derived from,
// so the gallery can never silently drift from the registered effect set when an
// effect is added or removed.
const REPO_ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');

async function loadEffectRoster() {
  const src = await readFile(join(REPO_ROOT, 'core', 'effects.h'), 'utf8');
  const block = src.match(/#define HS_EFFECT_LIST\(X\)([\s\S]*?)\r?\n\r?\n/);
  if (!block) throw new Error('Could not locate HS_EFFECT_LIST in core/effects.h');
  const names = [...block[1].matchAll(/X\((\w+)\)/g)].map(m => m[1]);
  if (names.length === 0) throw new Error('HS_EFFECT_LIST parsed to zero effects');
  return names;
}

const EFFECTS = await loadEffectRoster();

await mkdir(OUT_DIR, { recursive: true });

const browser = await chromium.launch({
  headless: true,
  args: [
    '--enable-webgl',
    '--use-gl=angle',
    '--use-angle=swiftshader',
    '--enable-unsafe-swiftshader',
    '--ignore-gpu-blocklist',
    '--enable-accelerated-2d-canvas',
  ],
});
const ctx = await browser.newContext({
  ignoreHTTPSErrors: true,
  viewport: { width: 1600, height: 1200 },
  deviceScaleFactor: 2,
});
// Force preserveDrawingBuffer:true on any WebGL context so we can read the
// drawn frame via canvas.toDataURL after rendering.
await ctx.addInitScript(() => {
  const origGet = HTMLCanvasElement.prototype.getContext;
  HTMLCanvasElement.prototype.getContext = function(type, attrs) {
    if (type === 'webgl' || type === 'webgl2' || type === 'experimental-webgl') {
      attrs = Object.assign({}, attrs || {}, { preserveDrawingBuffer: true });
    }
    return origGet.call(this, type, attrs);
  };
});

const page = await ctx.newPage();

page.on('console', msg => {
  const t = msg.type();
  if (t === 'error' || t === 'warning') console.log(`[${t}]`, msg.text());
});

const targets = process.argv.slice(2).length ? process.argv.slice(2) : EFFECTS;

let failures = 0;
for (const effect of targets) {
  const url = `${BASE_URL}?effect=${encodeURIComponent(effect)}&resolution=${encodeURIComponent(RESOLUTION)}`;
  process.stdout.write(`Capturing ${effect}... `);
  try {
    await page.goto(url, { waitUntil: 'load', timeout: 60000 });
    // Wait for canvas to exist and effect to settle/animate
    await page.waitForSelector('#canvas', { timeout: 30000 });
    await page.waitForTimeout(WAIT_MS);

    // With preserveDrawingBuffer:true forced via addInitScript, canvas
    // toDataURL is safe to call after rendering settles. Daydream's driver
    // suppresses the PiP under navigator.webdriver, so no post-crop needed.
    const dataUrl = await page.evaluate(() => {
      return document.querySelector('#canvas').toDataURL('image/png');
    });
    const b64 = dataUrl.split(',', 2)[1];
    const buf = Buffer.from(b64, 'base64');
    const out = `${OUT_DIR}/${effect}.png`;
    await writeFile(out, buf);
    console.log(`saved ${out}`);
  } catch (e) {
    failures++;
    console.log(`FAILED: ${e.message}`);
  }
}

await browser.close();

// A failed capture leaves the previous (stale) PNG in place, and that gallery is
// installed into daydream and served live, so a silent exit-0 would ship stale
// screenshots. Surface any failure to the caller (mirrors wasm_smoke.mjs).
if (failures) {
  console.log(`${failures} of ${targets.length} captures failed`);
  process.exitCode = 1;
}
