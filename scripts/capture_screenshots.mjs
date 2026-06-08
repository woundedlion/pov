// Headless Playwright script that loads the WebGL simulator for each effect,
// lets it animate, and saves a PNG screenshot to docs/screenshots/.
// Effects can be overridden via CLI args; otherwise the full EFFECTS list runs.
import { chromium } from 'playwright';
import { mkdir } from 'node:fs/promises';
import { dirname } from 'node:path';

const BASE_URL = process.env.SIM_URL || 'https://localhost/';
const RESOLUTION = 'Phantasm (144x288)';
const OUT_DIR = 'docs/screenshots';
const WAIT_MS = parseInt(process.env.WAIT_MS || '30000', 10);

const EFFECTS = [
  'BZReactionDiffusion', 'GSReactionDiffusion', 'HopfFibration', 'IslamicStars',
  'HankinSolids', 'SphericalHarmonics', 'Metaballs', 'MobiusGrid',
  'Moire', 'FlowField', 'Voronoi', 'PetalFlow', 'DreamBalls', 'Comets',
  'RingSpin', 'RingShower', 'ChaoticStrings', 'MeshFeedback', 'Liquid2D',
  'MindSplatter', 'Dynamo', 'Thrusters', 'GnomonicStars', 'Raymarch',
  'Flyby', 'SplineFlow',
];

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
    await (await import('node:fs/promises')).writeFile(out, buf);
    console.log(`saved ${out}`);
  } catch (e) {
    console.log(`FAILED: ${e.message}`);
  }
}

await browser.close();
