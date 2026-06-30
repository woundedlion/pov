// Headless Node smoke test for the shipped WASM module.
//
// The CI `wasm` job is otherwise compile-and-upload only: nothing instantiates
// the MODULARIZE/ES6 module, so a SIMD-codegen fault, an embind signature
// mismatch, a stack overflow (STACK_SIZE is a tight 8 KB), or an
// ALLOW_MEMORY_GROWTH detachment regression would pass CI and deploy. This
// script imports the built module exactly as a browser would, drives the engine
// across every registered effect at every enumerated resolution, and asserts the
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

const FRAMES_PER_EFFECT = 3;

// The stack is the one tracked region with no allocator trap: an overrun
// silently corrupts memory below stack-end, and stack_high_water_mark()
// saturates at `capacity` because the canary scan begins at stack-end and
// cannot see past it. So `hwm > capacity` (the arena check) is unsatisfiable
// for the stack — a fully consumed 8 KB stack reports hwm == capacity and would
// otherwise pass the exact gate CI claims catches a stack overflow. Require the
// stack to retain a real margin so a creeping regression trips long before the
// cliff. Current peak across all effects is ~0.6 KB, so 75% is ample headroom.
//
// Note this gate is a tripwire, not a guaranteed bound: high_water_mark is a
// LOWER bound that under-reports (a frame writing a coincidental canary byte, or
// reserving space it never stores to, reads back as still-canary — see
// stack_high_water_mark() in wasm.cpp). So a true stack usage above 75% can read
// back below the gate and pass; the 75% number catches creeping regressions, it
// does not certify the stack stayed under three-quarters full.
const STACK_MAX_FILL = 0.75;

// Wrapped in main() so a fatal precondition can `process.exitCode = 1; return`
// rather than process.exit(), letting buffered stdout/stderr flush before exit.
async function main() {
  try {
    await access(jsPath);
  } catch {
    console.error(`wasm_smoke: module not found at ${jsPath}\n` +
      `Build it first (cmake --preset wasm-release && cmake --build --preset wasm-release) ` +
      `or pass the path as an argument.`);
    process.exitCode = 1;
    return;
  }

  const { default: createHolosphereModule } = await import(pathToFileURL(jsPath));

  // Surface engine-side hs::log output and any abort() so a trap is visible in
  // the CI log rather than a bare non-zero exit.
  const Module = await createHolosphereModule({
    print: (s) => console.log(`[wasm] ${s}`),
    printErr: (s) => console.error(`[wasm:err] ${s}`),
  });

  let failures = 0;
  const fail = (msg) => { console.error(`  FAIL: ${msg}`); failures++; };

  // The supported (W,H) set is enumerated from the module itself (generated from
  // HS_WASM_RESOLUTIONS in wasm.cpp) rather than hand-mirrored here, so a newly
  // added resolution gets smoke coverage automatically — the same anti-drift
  // guarantee the effect roster already has via getEffectSizes().
  const RESOLUTIONS = Module.HolosphereEngine.getSupportedResolutions();
  if (!RESOLUTIONS || RESOLUTIONS.length === 0) {
    console.error('wasm_smoke: getSupportedResolutions() returned no resolutions');
    process.exitCode = 1;
    return;
  }

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
        // The stack traps nowhere, so guard it with a margin instead of waiting
        // for hwm > capacity (which it can never reach — see STACK_MAX_FILL).
        const stack = m.stack;
        if (!stack) {
          fail(`${name}: getArenaMetrics() omits the stack region`);
        } else if (stack.high_water_mark === 0) {
          fail(`${name}: stack high-water mark is 0 (canary tracking is broken)`);
        } else if (!(stack.capacity > 0)) {
          fail(`${name}: stack capacity ${stack.capacity} is not positive (the margin check below would divide a degenerate budget)`);
        } else if (stack.high_water_mark >= stack.capacity * STACK_MAX_FILL) {
          fail(`${name}: stack high-water mark ${stack.high_water_mark} of ` +
            `${stack.capacity} bytes leaves under ${Math.round((1 - STACK_MAX_FILL) * 100)}% ` +
            `margin — approaching overflow`);
        }

        // The embind param-marshalling layer — getParameterDefinitions() and the
        // per-frame getParamValues() stream — is the JS↔C++ order seam the GUI
        // rides every frame, and nothing else in this script exercises it. A
        // transposed binding, a wrong embind signature, or a length/order split
        // between the two would ride a green build straight to a desynced param
        // panel. Drive both here and assert they stay zippable and well-formed.
        const defs = engine.getParameterDefinitions();
        const values = engine.getParamValues();
        if (!Array.isArray(defs)) {
          fail(`${name}: getParameterDefinitions() did not return an array`);
        } else {
          // param_marshal.h is the single source of order shared by both calls.
          // A length split means the streams can no longer be zipped at all; the
          // per-index value check below then proves the surviving pairs are not
          // transposed — a reorder keeps lengths equal, so length alone is blind.
          if (values.length !== defs.length) {
            fail(`${name}: getParamValues() length ${values.length} != ` +
              `getParameterDefinitions() length ${defs.length} (param order seam drifted)`);
          }
          for (let i = 0; i < defs.length; i++) {
            const d = defs[i];
            if (typeof d.name !== 'string' || d.name.length === 0) {
              fail(`${name}: param ${i} has no name`);
            }
            // The value stream carries no names, so order is proven by checking
            // values[i] reproduces defs[i].value at the same index — both derive
            // from one Effect::getParameters() pass (param_marshal.h) with no
            // drawFrame between the two calls, so they are bit-identical per i.
            // A reordering bug keeps lengths equal but transposes the pairing,
            // tripping this where a length-only check stays green.
            const sv = i < values.length ? values[i] : undefined;
            if (i < values.length && !Number.isFinite(sv)) {
              fail(`${name}: param "${d.name}" value-stream entry ${sv} is not finite`);
            }
            if (typeof d.value === 'boolean') {
              // wasm.cpp collapses a bool def's value to `raw > 0.5`; the value
              // stream keeps the raw float, so reconstruct and compare.
              if (i < values.length && d.value !== (sv > 0.5)) {
                fail(`${name}: param "${d.name}" (index ${i}) def bool ${d.value} ` +
                  `!= value-stream ${sv} > 0.5 (param order seam transposed)`);
              }
              continue; // bools omit min/max (wasm.cpp)
            }
            // Float def: the stream entry must equal the def's value at this
            // index, or the two streams are not zippable.
            if (i < values.length && Number.isFinite(sv) && sv !== d.value) {
              fail(`${name}: param "${d.name}" (index ${i}) def value ${d.value} ` +
                `!= value-stream ${sv} (param order seam transposed)`);
            }
            // Float params carry a finite, ordered range bracketing their value.
            const eps = 1e-4 * (1 + Math.abs(d.max - d.min));
            if (!Number.isFinite(d.min) || !Number.isFinite(d.max) || d.min > d.max) {
              fail(`${name}: param "${d.name}" has a non-finite/inverted range [${d.min}, ${d.max}]`);
            } else if (!Number.isFinite(d.value) || d.value < d.min - eps || d.value > d.max + eps) {
              fail(`${name}: param "${d.name}" value ${d.value} outside [${d.min}, ${d.max}]`);
            }
          }
        }
      }
    }

    // Report the worst-case stack usage as a margin against STACK_SIZE so a
    // creeping regression is visible in the log before it ever overflows.
    const stack = engine.getArenaMetrics().stack;
    if (!stack) fail('getArenaMetrics() omits the stack region');
    else console.log(`\nstack: ${stack.high_water_mark}/${stack.capacity} bytes peak`);
  } finally {
    engine.delete();
  }

  // ── MeshOps tooling bindings + spline exports ───────────────────────────────
  // The engine loop above drives only HolosphereEngine. The MeshOpsWrapper surface
  // (lazy tooling arenas, the clearToolingMemory generation trap, the Conway
  // operators, getVertices/getFaces/classifyFaces) and the exported spline
  // functions are never instantiated — exactly the embind-signature drift this
  // script exists to catch, and a drift there rides a green build straight to the
  // live solids.html / splines.html tools. Exercise them here.
  console.log('\nMeshOps + splines:');

  const MeshOps = Module.MeshOps;
  if (!MeshOps) {
    fail('Module.MeshOps binding is missing');
  } else {
    // Use a real registry name (anti-drift, like the resolution/effect rosters
    // above) rather than hardcoding one.
    const registry = MeshOps.getRegistry();
    const solidName = registry && registry.length ? registry[0].name : null;
    if (!solidName) {
      fail('MeshOps.getRegistry() returned no solids');
    } else {
      // Unknown names must be rejected (null), not abort the module.
      const bogus = MeshOps.fromSolidName('definitely_not_a_solid');
      if (bogus) { fail('fromSolidName(unknown) should return null'); bogus.delete(); }

      const solid = MeshOps.fromSolidName(solidName);
      if (!solid) {
        fail(`fromSolidName("${solidName}") returned null`);
      } else {
        const verts = solid.getVertices();
        if (!verts || verts.length === 0 || verts.length % 3 !== 0) {
          fail(`${solidName} getVertices(): ${verts && verts.length} floats (want nonzero multiple of 3)`);
        }
        if (solid.getFaces().length === 0) fail(`${solidName} getFaces() returned no faces`);

        // Conway operator path: dual() runs apply()/finalize into the tooling arena.
        const dual = solid.dual();
        if (!dual) {
          fail(`${solidName}.dual() returned null`);
        } else {
          if (dual.getVertices().length === 0) fail('dual getVertices() empty');
          if (dual.classifyFaces().length === 0) fail('dual classifyFaces() empty');
          dual.delete();
        }

        // The parameterized embind operators — the float-arg MESHOP_1F path, the
        // int-arg relax with its kMaxRelaxIterations clamp, and the finite-arg
        // hankin reject — are exactly the arg-marshaling/signature seams this
        // script exists to catch, and none are exercised above. Drive a
        // representative of each on the live solid.
        const isValidMesh = (w) =>
          w && w.getVertices().length > 0 && w.getVertices().length % 3 === 0 &&
          w.getFaces().length > 0;

        // MESHOP_1F(truncate): marshals a float arg and finalizes a new mesh.
        const trunc = solid.truncate(0.3);
        if (!isValidMesh(trunc)) fail(`${solidName}.truncate(0.3) did not produce a valid mesh`);
        if (trunc) trunc.delete();

        // relax(int): the int-arg path plus its C++-side clamp. relax(1) is the
        // nominal case; relax(1e9) is a large-but-still-INT32-valid count
        // (1e9 < 2^31) that must clamp to kMaxRelaxIterations rather than loop a
        // billion times, and still return a valid mesh. This exercises the C++
        // clamp only — NOT embind's double->int coercion at/past INT32_MAX, which
        // is a separate marshaling path this test does not reach.
        const relaxed = solid.relax(1);
        if (!isValidMesh(relaxed)) fail(`${solidName}.relax(1) did not produce a valid mesh`);
        if (relaxed) relaxed.delete();
        const relaxedCap = solid.relax(1e9);
        if (!isValidMesh(relaxedCap)) fail(`${solidName}.relax(1e9) did not clamp to a valid mesh`);
        if (relaxedCap) relaxedCap.delete();

        // hankin(float): a non-finite arg must be rejected at the boundary
        // (finite_arg → null) rather than abort the module.
        const hankinBad = solid.hankin(NaN);
        if (hankinBad) { fail(`${solidName}.hankin(NaN) should return null`); hankinBad.delete(); }

        solid.delete();
      }

      // Tooling arena high-water marks must stay within capacity after the ops.
      const tm = MeshOps.getArenaMetrics();
      for (const region of Object.keys(tm)) {
        const { high_water_mark: hwm, capacity } = tm[region];
        if (hwm > capacity) fail(`MeshOps ${region} high-water mark ${hwm} exceeds capacity ${capacity}`);
      }

      // The clearToolingMemory generation trap: the wipe reclaims the tooling
      // arenas; a fresh build afterwards must still succeed.
      MeshOps.clearToolingMemory();
      const post = MeshOps.fromSolidName(solidName);
      if (!post) {
        fail('fromSolidName after clearToolingMemory() returned null');
      } else {
        if (post.getVertices().length === 0) fail('post-wipe getVertices() empty');
        post.delete();
      }
      console.log(`  MeshOps: ${solidName} + dual, truncate, relax(+clamp), hankin reject, classifyFaces, clearToolingMemory OK`);
    }
  }

  // Spline exports used by splines.html. Assert finite {x,y,z} AND pin numeric
  // behavior, so a transposed-argument or wrong-target embind binding fails — not
  // only a non-finite return.
  const isVec = (v) => v && Number.isFinite(v.x) && Number.isFinite(v.y) && Number.isFinite(v.z);
  const approxVec = (v, x, y, z, eps = 1e-4) =>
    isVec(v) && Math.abs(v.x - x) <= eps && Math.abs(v.y - y) <= eps && Math.abs(v.z - z) <= eps;

  // p0..p3 control polygon (all unit vectors), shared by the cubic evaluators.
  const CTRL = [0, 1, 0, 1, 0, 0, 0, 0, 1, -1, 0, 0];
  for (const name of ['spline_cubic_fast', 'spline_cubic_slerp']) {
    const mid = Module[name](...CTRL, 0.5);
    if (!isVec(mid)) fail(`${name} returned non-finite ${JSON.stringify(mid)}`);
    // A cubic evaluator must interpolate its endpoints: p0 at t=0, p3 at t=1.
    // Both control points are unit, so the normalized endpoints are p0/p3 exactly.
    // Catches argument transposition and a binding wired to the wrong output.
    const at0 = Module[name](...CTRL, 0);
    const at1 = Module[name](...CTRL, 1);
    if (!approxVec(at0, 0, 1, 0)) fail(`${name}(t=0) should be p0=(0,1,0), got ${JSON.stringify(at0)}`);
    if (!approxVec(at1, -1, 0, 0)) fail(`${name}(t=1) should be p3=(-1,0,0), got ${JSON.stringify(at1)}`);
  }
  // cubic_fast (normalized Bézier) and cubic_slerp (spherical de Casteljau) are
  // genuinely different curves away from the endpoints; identical midpoints would
  // mean both bindings resolve to the same underlying function (embind drift).
  const fastMid = Module.spline_cubic_fast(...CTRL, 0.5);
  const slerpMid = Module.spline_cubic_slerp(...CTRL, 0.5);
  if (approxVec(slerpMid, fastMid.x, fastMid.y, fastMid.z)) {
    fail(`cubic_fast and cubic_slerp produced identical midpoints ${JSON.stringify(fastMid)} — bindings may resolve to the same function`);
  }
  // prev=(0,1,0) start=(1,0,0) end=(0,0,1) next=(-1,0,0), tension=0.5.
  const tg = Module.spline_catmull_rom_tangents(0, 1, 0, 1, 0, 0, 0, 0, 1, -1, 0, 0, 0.5);
  if (!tg || !isVec(tg.cp1) || !isVec(tg.cp2)) {
    fail(`spline_catmull_rom_tangents returned malformed ${JSON.stringify(tg)}`);
  }
  // Pin cp1 so an argument transposition or a wrong-output binding fails here, not
  // only a non-finite return. cp1 = slerp(start, slerp(prev,end,.5), .5); both
  // slerps join orthonormal vectors at t=0.5, where the equal weights cancel under
  // the final normalize, so cp1 = normalize(start + normalize(prev+end)) =
  // (√½, ½, ½) exactly — independent of the fast-trig approximation.
  else if (!approxVec(tg.cp1, Math.SQRT1_2, 0.5, 0.5)) {
    fail(`spline_catmull_rom_tangents cp1 should be (0.7071, 0.5, 0.5), got ${JSON.stringify(tg.cp1)}`);
  }
  console.log('  splines: cubic_fast, cubic_slerp, catmull_rom_tangents OK');

  if (failures > 0) {
    console.error(`\nwasm_smoke: ${failures} failure(s)`);
    process.exitCode = 1;
    return;
  }
  console.log('\nwasm_smoke: OK');
}

await main();
