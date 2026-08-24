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
import { mkdir, writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import { loadEffectRoster, REPO_ROOT } from './effect_roster.mjs';
import {
  captureIsBlank,
  captureLitFraction,
  captureOffsetMs,
  COVERAGE_SAMPLE_HEIGHT,
  COVERAGE_SAMPLE_WIDTH,
  DEFAULT_BLANK_FLOOR,
  DEFAULT_CAPTURE_OFFSET_MS,
  GALLERY_HEIGHT,
  GALLERY_WIDTH,
  minimumLitPixels,
} from './screenshot_capture_config.mjs';
import { descendToHonoredResolution } from './screenshot_resolution.mjs';

// Number('') is 0 (finite), so blank/whitespace is rejected explicitly.
async function numEnv(name, def) {
  const raw = process.env[name];
  if (raw === undefined) return def;
  const v = Number(raw);
  if (raw.trim() !== '' && Number.isFinite(v) && v >= 0) return v;
  console.error('========================================================');
  console.error(`capture_screenshots: ERROR — ${name} must be a finite, non-negative number.`);
  console.error(`Received: ${JSON.stringify(raw)}`);
  console.error('========================================================');
  process.exitCode = 2;
  // Drain buffered stderr before the hard exit; a pipe truncates it otherwise.
  await new Promise((resolve) => process.stderr.write('', resolve));
  process.exit();
}

const BASE_URL = process.env.SIM_URL || 'http://localhost:8080/';
const OUT_DIR = join(REPO_ROOT, 'docs', 'screenshots');
const WAIT_MS = await numEnv('WAIT_MS', DEFAULT_CAPTURE_OFFSET_MS);
const WAIT_MS_OVERRIDE = process.env.WAIT_MS === undefined ? null : WAIT_MS;
const BLANK_FLOOR = await numEnv('BLANK_FLOOR', DEFAULT_BLANK_FLOOR);

const { chromium } = await import('playwright');

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
// Thrown when the page yields no resolution list, to abort the capture run from
// inside the browser block without leaking a raw stack past the summaries.
class UnresolvedResolutions extends Error {}

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
  // filename. On any failure, return [] — the caller aborts the run rather than
  // capturing unpinned.
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
      console.warn('capture_screenshots: found no resolution options on the page');
    } catch (e) {
      console.warn('capture_screenshots: could not read resolutions from the page ' +
        `(${e.message})`);
    }
    return [];
  }

  // Resolutions to try per effect, largest first: each is tried from highest to
  // lowest detail, keeping the first that honors the requested effect.
  RESOLUTIONS = await resolveResolutions();
  // Without a resolution list every URL would omit the resolution param, and a
  // request with no param cannot be confirmed against the app's rewritten effect
  // param — the app's fallback effect would be saved under the requested effect's
  // filename. Abort here, before any PNG is written.
  if (RESOLUTIONS.length === 0) throw new UnresolvedResolutions();
  console.log(`Capture resolutions (high→low): ${RESOLUTIONS.join(', ')}`);

  targets = process.argv.slice(2).length ? process.argv.slice(2) : EFFECTS;

  // Grab the current #canvas frame and measure how much of it is lit. With
  // preserveDrawingBuffer:true forced via addInitScript, toDataURL is safe after
  // rendering settles. Coverage is measured on a small downscale (cheap, and the
  // thumbnail is downscaled anyway): the fraction of pixels above a near-black
  // floor. Daydream's driver suppresses the PiP under navigator.webdriver, so no
  // post-crop is needed.
  async function grabFrame() {
    return await page.evaluate(({
      outWidth,
      outHeight,
      sampleWidth,
      sampleHeight,
    }) => {
      const canvas = document.querySelector('#canvas');
      const off = document.createElement('canvas');
      off.width = sampleWidth; off.height = sampleHeight;
      const g = off.getContext('2d');
      g.drawImage(canvas, 0, 0, sampleWidth, sampleHeight);
      const data = g.getImageData(0, 0, sampleWidth, sampleHeight).data;
      let lit = 0;
      for (let i = 0; i < data.length; i += 4) {
        if (data[i] > 12 || data[i + 1] > 12 || data[i + 2] > 12) lit++;
      }
      const thumb = document.createElement('canvas');
      thumb.width = outWidth;
      thumb.height = outHeight;
      thumb.getContext('2d').drawImage(canvas, 0, 0, thumb.width, thumb.height);
      return { dataUrl: thumb.toDataURL('image/png'), litPixels: lit };
    }, {
      outWidth: GALLERY_WIDTH,
      outHeight: GALLERY_HEIGHT,
      sampleWidth: COVERAGE_SAMPLE_WIDTH,
      sampleHeight: COVERAGE_SAMPLE_HEIGHT,
    });
  }

  // Loads one effect at one resolution and reports the effect the app actually
  // selected: it rewrites the URL's effect param to its choice, so a silent
  // fallback (requested effect not offered here) shows up as a different name.
  // descendToHonoredResolution() owns which resolution wins.
  async function loadEffect(effect, resolution) {
    const params = new URLSearchParams({ effect, resolution });
    await page.goto(`${BASE_URL}?${params.toString()}`,
      { waitUntil: 'load', timeout: 60000 });
    await page.waitForSelector('#canvas', { timeout: 30000 });
    // The fallback rewrite happens during hydration, before the settle wait.
    await page.waitForTimeout(500);
    return await page.evaluate(() =>
      new URLSearchParams(location.search).get('effect'));
  }

  for (const effect of targets) {
    process.stdout.write(`Capturing ${effect}... `);
    try {
      // Try resolutions high→low; keep the first that actually offers this effect.
      const { resolution: usedRes, honored } =
        await descendToHonoredResolution(effect, RESOLUTIONS, loadEffect);
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

      const { dataUrl, litPixels } = await grabFrame();
      const lit = captureLitFraction(litPixels);
      const pct = (lit * 100).toFixed(2);
      if (captureIsBlank(litPixels, BLANK_FLOOR)) {
        blanks.push(effect);
        console.log(
          `SKIPPED — fixed ${offsetMs}ms capture was blank ` +
            `(${litPixels} lit pixels; minimum ${minimumLitPixels(BLANK_FLOOR)}); ` +
            'kept existing PNG',
        );
        continue;
      }

      const b64 = dataUrl.split(',', 2)[1];
      const buf = Buffer.from(b64, 'base64');
      const out = join(OUT_DIR, `${effect}.png`);
      await writeFile(out, buf);
      console.log(`saved ${out} @ ${usedRes} after ${offsetMs}ms ` +
        `(${pct}% lit)`);
    } catch (e) {
      failures++;
      console.log(`FAILED: ${e.message}`);
    }
  }

} catch (e) {
  if (!(e instanceof UnresolvedResolutions)) throw e;
} finally {
  if (browser) await browser.close();
}

// resolveResolutions() returned [], so the run aborted before capturing: with no
// resolution param the app's silent fallback cannot be detected, and every PNG
// would risk carrying the fallback effect under another effect's filename. Its
// per-failure console.warn fires up front and scrolls away, so restate it in the
// summary the caller actually reads.
if (RESOLUTIONS.length === 0) {
  console.warn('========================================================');
  console.warn('capture_screenshots: ERROR — resolutions were NOT resolved;');
  console.warn('ABORTED before capturing. No screenshot was written; the existing');
  console.warn('gallery is untouched. Re-run once the page resolution selector is');
  console.warn('reachable.');
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
