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
// http.server port); WAIT_MS overrides every configured capture offset.
import { chromium } from 'playwright';
import { mkdir, writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import { loadEffectRoster, REPO_ROOT } from './effect_roster.mjs';
import { captureOffsetMs } from './screenshot_capture_config.mjs';

// A malformed env value silently disables or distorts the timing it controls,
// so fall back to the default on anything that isn't a finite, non-negative
// number. Number('') is 0 (finite), so blank/whitespace is rejected explicitly.
function numEnv(name, def) {
  const raw = process.env[name];
  if (raw === undefined || raw.trim() === '') return def;
  const v = Number(raw);
  return Number.isFinite(v) && v >= 0 ? v : def;
}

const BASE_URL = process.env.SIM_URL || 'http://localhost:8080/';
const OUT_DIR = join(REPO_ROOT, 'docs', 'screenshots');
const WAIT_MS_OVERRIDE = process.env.WAIT_MS === undefined
  ? null
  : numEnv('WAIT_MS', null);
const BLANK_FLOOR = numEnv('BLANK_FLOOR', 0.0005);

// The effect roster (and the docs/screenshots freshness gate that mirrors it)
// is parsed from the HS_EFFECT_LIST X-macro by scripts/effect_roster.mjs.
const EFFECTS = await loadEffectRoster();

await mkdir(OUT_DIR, { recursive: true });

// A launch failure (browser not installed) would otherwise throw a raw stack
// past the script's banner summaries; report it through the same actionable path.
let browser;
try {
  browser = await chromium.launch({
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
} catch (e) {
  console.warn('========================================================');
  console.warn(`capture_screenshots: ERROR — could not launch Chromium (${e.message}).`);
  console.warn('Install the browser once with:  npx playwright install chromium');
  console.warn('========================================================');
  process.exitCode = 1;
  // Drain buffered stderr before the hard exit; a pipe truncates it otherwise.
  await new Promise((resolve) => process.stderr.write('', resolve));
  process.exit();
}
// Declared out here, not inside the try below: the summary/gating section after
// the finally reads them, so block-scoping them to the try would leave every run
// throwing a ReferenceError past browser.close().
let RESOLUTIONS = [];
let targets = [];
let failures = 0;
const blanks = [];
const wrongRes = [];
try {
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

  // Read the app's supported resolutions instead of hard-coding a display label
  // (the same anti-drift goal as the effect roster above). The resolution
  // control's options are built from the app's supported resolutions; return them
  // sorted by pixel area (largest/highest-detail first). Each effect is then
  // captured at the FIRST resolution that actually offers it — not blindly at the
  // global largest, because the effect set is resolution-specific (e.g. RingShower
  // and Dynamo are only registered at the low-res Holosphere preset). Requesting
  // an effect at a resolution that doesn't offer it makes the app silently fall
  // back to its default effect, which would save that default under the wrong
  // filename. On any failure, return [] so each effect is captured with no
  // resolution param (the app picks its own default).
  async function resolveResolutions() {
    try {
      await page.goto(BASE_URL, { waitUntil: 'load', timeout: 60000 });
      await page.waitForSelector('select', { timeout: 30000 });
      const labels = await page.evaluate(() => {
        const re = /\((\d+)x(\d+)\)/; // resolution labels embed their dimensions
        const seen = new Map();
        for (const opt of document.querySelectorAll('option')) {
          const text = (opt.textContent || '').trim();
          const m = re.exec(text);
          if (m) seen.set(text, Number(m[1]) * Number(m[2]));
        }
        return [...seen.entries()].sort((a, b) => b[1] - a[1]).map(e => e[0]);
      });
      if (labels.length) return labels;
      console.warn('capture_screenshots: found no resolution options on the page; ' +
        'using the app default resolution');
    } catch (e) {
      console.warn(`capture_screenshots: could not read resolutions from the page ` +
        `(${e.message}); using the app default resolution`);
    }
    return [];
  }

  // Resolutions to try per effect, largest first; [null] means "no param, app
  // default". When the app exposes its presets we try each from highest to lowest
  // detail and keep the first that honors the requested effect.
  RESOLUTIONS = await resolveResolutions();
  const RES_TRY = RESOLUTIONS.length ? RESOLUTIONS : [null];
  console.log(`Capture resolutions (high→low): ${RESOLUTIONS.join(', ') || '(app default)'}`);

  targets = process.argv.slice(2).length ? process.argv.slice(2) : EFFECTS;

  // Grab the current #canvas frame and measure how much of it is lit. With
  // preserveDrawingBuffer:true forced via addInitScript, toDataURL is safe after
  // rendering settles. Coverage is measured on a small downscale (cheap, and the
  // thumbnail is downscaled anyway): the fraction of pixels above a near-black
  // floor. Daydream's driver suppresses the PiP under navigator.webdriver, so no
  // post-crop is needed.
  async function grabFrame() {
    return await page.evaluate(() => {
      const canvas = document.querySelector('#canvas');
      const SW = 96, SH = 72;
      const off = document.createElement('canvas');
      off.width = SW; off.height = SH;
      const g = off.getContext('2d');
      g.drawImage(canvas, 0, 0, SW, SH);
      const data = g.getImageData(0, 0, SW, SH).data;
      let lit = 0;
      for (let i = 0; i < data.length; i += 4) {
        if (data[i] > 12 || data[i + 1] > 12 || data[i + 2] > 12) lit++;
      }
      return { dataUrl: canvas.toDataURL('image/png'), lit: lit / (SW * SH) };
    });
  }

  // The app rewrites the URL's effect param to whatever it actually selected, so
  // after navigating we can detect a silent fallback (requested effect not offered
  // at this resolution) by comparing the rewritten param to what we asked for.
  async function selectedEffect() {
    return await page.evaluate(() =>
      new URLSearchParams(location.search).get('effect'));
  }

  for (const effect of targets) {
    process.stdout.write(`Capturing ${effect}... `);
    try {
      // Try resolutions high→low; keep the first that actually offers this effect.
      let usedRes = null, honored = false;
      for (const res of RES_TRY) {
        const params = new URLSearchParams({ effect });
        if (res) params.set('resolution', res);
        await page.goto(`${BASE_URL}?${params.toString()}`,
          { waitUntil: 'load', timeout: 60000 });
        await page.waitForSelector('#canvas', { timeout: 30000 });
        // The fallback rewrite happens during hydration, before the settle wait.
        await page.waitForTimeout(500);
        usedRes = res;
        if (res === null || (await selectedEffect()) === effect) { honored = true; break; }
      }
      // Offered at no resolution: the canvas shows the app's fallback effect.
      // Saving it would overwrite a (possibly correct) existing PNG with a
      // thumbnail of the WRONG effect — worse than leaving the stale one. Skip the
      // save and flag it; the prior PNG stays untouched.
      if (!honored) {
        wrongRes.push(effect);
        console.log(`SKIPPED — offered at no resolution (app fell back); kept existing PNG`);
        continue;
      }

      const offsetMs = captureOffsetMs(effect, WAIT_MS_OVERRIDE);
      await page.waitForTimeout(offsetMs);

      const { dataUrl, lit } = await grabFrame();
      const pct = (lit * 100).toFixed(2);
      if (lit < BLANK_FLOOR) {
        blanks.push(effect);
        console.log(`SKIPPED — fixed ${offsetMs}ms capture was blank (${pct}% lit); ` +
          'kept existing PNG');
        continue;
      }

      const b64 = dataUrl.split(',', 2)[1];
      const buf = Buffer.from(b64, 'base64');
      const out = join(OUT_DIR, `${effect}.png`);
      await writeFile(out, buf);
      console.log(`saved ${out} @ ${usedRes || 'default'} after ${offsetMs}ms ` +
        `(${pct}% lit)`);
    } catch (e) {
      failures++;
      console.log(`FAILED: ${e.message}`);
    }
  }

} finally {
  if (browser) await browser.close();
}

// resolveResolutions() returned [], so the per-capture URLs omitted the
// resolution param and the WHOLE gallery was captured at whatever the app
// defaulted to — not a pinned resolution. resolveResolutions deliberately
// degrades rather than aborting the run, but its per-failure console.warn fires
// up front and scrolls away in a long capture, so restate it loudly in the
// summary the caller actually reads: a gallery shipped at the wrong resolution
// must not look like a clean success.
if (RESOLUTIONS.length === 0) {
  console.warn('========================================================');
  console.warn('capture_screenshots: WARNING — resolutions were NOT resolved;');
  console.warn('the entire gallery was captured at the APP DEFAULT resolution');
  console.warn('(unverified). Re-run once the page resolution selector is');
  console.warn('reachable to capture at a pinned resolution.');
  console.warn('========================================================');
  process.exitCode = 1;
}

// An effect that the app offered at no available resolution was SKIPPED (its
// existing PNG left untouched) rather than overwritten with the fallback effect.
// The effect is registered in the roster but absent from the app's per-resolution
// effect lists, so the live app cannot select it either — surface it loudly and
// fail the run.
if (wrongRes.length) {
  console.warn('========================================================');
  console.warn(`capture_screenshots: WARNING — ${wrongRes.length} effect(s) offered at NO`);
  console.warn('available resolution; SKIPPED (kept existing PNG, not regenerated):');
  console.warn(`  ${wrongRes.join(', ')}`);
  console.warn('These are in the roster but absent from the app\'s per-resolution');
  console.warn('effect lists, so the live app cannot select them. Add them to a');
  console.warn('resolution\'s effect list (daydream) or remove them from the roster.');
  console.warn('========================================================');
  process.exitCode = 1;
}

// A blank fixed-offset frame is not saved, so the previous capture remains.
if (blanks.length) {
  console.warn('========================================================');
  console.warn(`capture_screenshots: WARNING — ${blanks.length} capture(s) were BLANK:`);
  console.warn(`  ${blanks.join(', ')}`);
  console.warn('Adjust its fixed capture offset, or check the effect.');
  console.warn('========================================================');
  process.exitCode = 1;
}

// A failed capture leaves the previous (stale) PNG in place, and that gallery is
// installed into daydream and served live, so a silent exit-0 would ship stale
// screenshots. Surface any failure to the caller (mirrors wasm_smoke.mjs).
if (failures) {
  console.log(`${failures} of ${targets.length} captures failed`);
  process.exitCode = 1;
}
