// Headless Node smoke test for the shipped WASM module.
//
// The CI `wasm` job is otherwise compile-and-upload only: nothing instantiates
// the MODULARIZE/ES6 module, so a SIMD-codegen fault, an embind signature
// mismatch, a stack overflow (STACK_SIZE is a tight 8 KB), or an
// ALLOW_MEMORY_GROWTH detachment regression would pass CI and deploy. This
// script imports the built module exactly as a browser would, drives the engine
// across every registered effect at both supported resolutions, and asserts the
// arena/stack high-water marks stay within their capacities — turning a runtime
// regression into a red build instead of a broken deploy.
//
// Usage (from the Holosphere repo root, after a wasm build):
//   node scripts/wasm_smoke.mjs [path/to/holosphere_wasm.js]
// Defaults to the wasm-release build output; override with the arg or WASM_JS.
import { pathToFileURL } from 'node:url';
import { join, isAbsolute } from 'node:path';
import { access } from 'node:fs/promises';

const DEFAULT_JS = 'build/wasm-release/holosphere_wasm.js';
const jsArg = process.argv[2] || process.env.WASM_JS || DEFAULT_JS;
const jsPath = isAbsolute(jsArg) ? jsArg : join(process.cwd(), jsArg);

// Resolutions the WASM factory supports (mirrors HS_WASM_RESOLUTIONS in
// targets/wasm/wasm.cpp). Each effect is exercised at both.
const RESOLUTIONS = [
  [96, 20],
  [288, 144],
];
const FRAMES_PER_EFFECT = 3;

await access(jsPath).catch(() => {
  console.error(`wasm_smoke: module not found at ${jsPath}\n` +
    `Build it first (cmake --preset wasm-release && cmake --build --preset wasm-release) ` +
    `or pass the path as an argument.`);
  process.exit(1);
});

const { default: createHolosphereModule } = await import(pathToFileURL(jsPath));

// Surface engine-side hs::log output and any abort() so a trap is visible in
// the CI log rather than a bare non-zero exit.
const Module = await createHolosphereModule({
  print: (s) => console.log(`[wasm] ${s}`),
  printErr: (s) => console.error(`[wasm:err] ${s}`),
});

let failures = 0;
const fail = (msg) => { console.error(`  FAIL: ${msg}`); failures++; };

const engine = new Module.HolosphereEngine();
try {
  for (const [w, h] of RESOLUTIONS) {
    if (!engine.setResolution(w, h)) {
      fail(`setResolution(${w}, ${h}) rejected a supported resolution`);
      continue;
    }

    // The running registry is the single source of truth for the effect roster;
    // enumerating it here means the smoke test can never drift from the set the
    // module actually built (same guarantee the screenshot gallery relies on).
    const sizes = engine.getEffectSizes();
    const names = Object.keys(sizes);
    if (names.length === 0) {
      fail(`getEffectSizes() returned no effects at ${w}x${h}`);
      continue;
    }
    console.log(`\n${w}x${h}: ${names.length} effects`);

    for (const name of names) {
      if (!engine.setEffect(name)) {
        fail(`setEffect("${name}") returned false at ${w}x${h}`);
        continue;
      }
      for (let f = 0; f < FRAMES_PER_EFFECT; f++) engine.drawFrame();

      // getPixels() aliases WASM memory; after drawing it must expose the full
      // active-resolution RGB span and not be a detached zero-length view.
      const px = engine.getPixels();
      const expected = w * h * 3;
      if (px.length !== expected) {
        fail(`${name}: getPixels() length ${px.length}, expected ${expected} ` +
          `(detached view or wrong stride)`);
      }

      // The whole point of the test: prove the tight 8 KB stack and the arenas
      // were not overrun while rendering this effect. The module reports each
      // region's high-water mark and capacity, so compare them directly — no
      // need to hardcode STACK_SIZE here (the build owns that number).
      const m = engine.getArenaMetrics();
      for (const region of Object.keys(m)) {
        const { high_water_mark: hwm, capacity } = m[region];
        if (hwm > capacity) {
          fail(`${name}: ${region} high-water mark ${hwm} exceeds capacity ${capacity}`);
        }
      }
      if (m.stack.high_water_mark === 0) {
        fail(`${name}: stack high-water mark is 0 (canary tracking is broken)`);
      }
    }
  }

  // Report the worst-case stack usage as a margin against STACK_SIZE so a
  // creeping regression is visible in the log before it ever overflows.
  const stack = engine.getArenaMetrics().stack;
  console.log(`\nstack: ${stack.high_water_mark}/${stack.capacity} bytes peak`);
} finally {
  engine.delete();
}

if (failures > 0) {
  console.error(`\nwasm_smoke: ${failures} failure(s)`);
  process.exit(1);
}
console.log('\nwasm_smoke: OK');
