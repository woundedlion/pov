// Headless Node smoke test for the shipped WASM module.
//
// The CI `wasm` job is otherwise compile-and-upload only, so a SIMD/codegen
// fault, embind signature mismatch, stack overflow, or memory-growth detachment
// would deploy unseen. Drives every effect at every resolution and asserts the
// arena/stack high-water marks stay within capacity.
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

// CI overrides via WASM_SMOKE_FRAMES to reach late-lifecycle events the
// 3-frame default never hits (frame-48 ShapeShifter cut, arena compaction).
const FRAMES_PER_EFFECT = Number(process.env.WASM_SMOKE_FRAMES ?? 3);

// The stack has no allocator trap and stack_high_water_mark() saturates at
// capacity, so `hwm > capacity` can never fire for it. Gate on an absolute byte
// budget instead (meaningful against any build's stack size); the min with the
// capacity fraction covers a hypothetical sub-ceiling stack. Creep tripwire, not
// a bound: the HWM under-reports (see stack_high_water_mark() in wasm.cpp).
// The 2048 default is calibrated on the -O3 release build; -O0 debug frames run
// severalfold larger, so ci.yml overrides via WASM_SMOKE_STACK_CEILING.
const STACK_HWM_CEILING_BYTES = Number(process.env.WASM_SMOKE_STACK_CEILING ?? 2048);
const STACK_MAX_FILL = 0.75;

// main() lets a fatal precondition set exitCode and return, so buffered stdout
// flushes rather than being cut off by process.exit().
async function main() {
  if (!Number.isInteger(STACK_HWM_CEILING_BYTES) || STACK_HWM_CEILING_BYTES <= 0) {
    console.error(`wasm_smoke: WASM_SMOKE_STACK_CEILING must be a positive integer, ` +
      `got "${process.env.WASM_SMOKE_STACK_CEILING}"`);
    process.exitCode = 1;
    return;
  }
  if (!Number.isInteger(FRAMES_PER_EFFECT) || FRAMES_PER_EFFECT <= 0) {
    console.error(`wasm_smoke: WASM_SMOKE_FRAMES must be a positive integer, ` +
      `got "${process.env.WASM_SMOKE_FRAMES}"`);
    process.exitCode = 1;
    return;
  }
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

  // Enumerated from the module (generated from HS_WASM_RESOLUTIONS) so a new
  // resolution gets coverage without editing this file.
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

      // Enumerate from the running registry so the roster can't drift from the
      // set the module actually built.
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

        // Assert no arena was overrun rendering this effect; the module reports
        // each region's high-water mark and capacity.
        const m = engine.getArenaMetrics();
        for (const region of Object.keys(m)) {
          const { high_water_mark: hwm, capacity } = m[region];
          if (hwm > capacity) {
            fail(`${name}: ${region} high-water mark ${hwm} exceeds capacity ${capacity}`);
          }
        }
        // The stack traps nowhere: guard it with the creep budget, not
        // hwm > capacity (unreachable — see STACK_HWM_CEILING_BYTES).
        const stack = m.stack;
        const stackGate = Math.min(STACK_HWM_CEILING_BYTES, stack && stack.capacity > 0
          ? stack.capacity * STACK_MAX_FILL : STACK_HWM_CEILING_BYTES);
        if (!stack) {
          fail(`${name}: getArenaMetrics() omits the stack region`);
        } else if (stack.high_water_mark === 0) {
          fail(`${name}: stack high-water mark is 0 (canary tracking is broken)`);
        } else if (!(stack.capacity > 0)) {
          fail(`${name}: stack capacity ${stack.capacity} is not positive (the gate below would compare against a degenerate budget)`);
        } else if (stack.high_water_mark >= stackGate) {
          fail(`${name}: stack high-water mark ${stack.high_water_mark} of ` +
            `${stack.capacity} bytes exceeds the ${stackGate}-byte creep budget — approaching overflow`);
        }

        // Exercise the embind param seam (getParameterDefinitions() +
        // getParamValues()) the GUI rides every frame: assert the two streams
        // stay zippable and well-formed.
        const defs = engine.getParameterDefinitions();
        const values = engine.getParamValues();
        if (!Array.isArray(defs)) {
          fail(`${name}: getParameterDefinitions() did not return an array`);
        } else {
          // Both calls share param_marshal.h's ordering. Length guards a split;
          // the per-index check below catches a reorder (length alone is blind
          // to transposition).
          if (values.length !== defs.length) {
            fail(`${name}: getParamValues() length ${values.length} != ` +
              `getParameterDefinitions() length ${defs.length} (param order seam drifted)`);
          }
          for (let i = 0; i < defs.length; i++) {
            const d = defs[i];
            if (typeof d.name !== 'string' || d.name.length === 0) {
              fail(`${name}: param ${i} has no name`);
            }
            // The value stream carries no names; values[i] must reproduce
            // defs[i].value (same getParameters() pass, no drawFrame between),
            // so a transposition trips here.
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

    // Log worst-case stack usage as a margin against STACK_SIZE.
    const stack = engine.getArenaMetrics().stack;
    if (!stack) fail('getArenaMetrics() omits the stack region');
    else console.log(`\nstack: ${stack.high_water_mark}/${stack.capacity} bytes peak`);
  } finally {
    engine.delete();
  }

  // ── MeshOps tooling bindings + spline exports ───────────────────────────────
  // The engine loop drives only HolosphereEngine; exercise the MeshOpsWrapper
  // surface and exported spline functions (used by solids.html / splines.html)
  // so their embind signatures can't drift unseen.
  console.log('\nMeshOps + splines:');

  const MeshOps = Module.MeshOps;
  if (!MeshOps) {
    fail('Module.MeshOps binding is missing');
  } else {
    // Use a real registry name rather than hardcoding one (anti-drift).
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
        // getFaces() returns parallel flat buffers: { indices, counts }.
        const faces = solid.getFaces();
        if (!faces || !faces.counts || faces.counts.length === 0) {
          fail(`${solidName} getFaces() returned no faces`);
        } else if (faces.indices.length !== faces.counts.reduce((sum, n) => sum + n, 0)) {
          fail(`${solidName} getFaces(): indices length ${faces.indices.length} ` +
            `disagrees with the counts sum`);
        }

        // Conway operator path: dual() runs apply()/finalize into the tooling arena.
        const dual = solid.dual();
        if (!dual) {
          fail(`${solidName}.dual() returned null`);
        } else {
          if (dual.getVertices().length === 0) fail('dual getVertices() empty');
          if (dual.classifyFaces().length === 0) fail('dual classifyFaces() empty');
          dual.delete();
        }

        // Drive the parameterized operators (float-arg truncate, int-arg relax
        // with its clamp, finite-arg hankin reject) — arg-marshaling seams not
        // exercised above.
        const isValidMesh = (w) => {
          if (!w) return false;
          const v = w.getVertices();
          const f = w.getFaces();
          return v.length > 0 && v.length % 3 === 0 &&
            f && f.counts && f.counts.length > 0 &&
            f.indices.length === f.counts.reduce((sum, n) => sum + n, 0);
        };

        // MESHOP_1U(truncate): marshals a [0,1]-clamped arg and finalizes a new mesh.
        const trunc = solid.truncate(0.3);
        if (!isValidMesh(trunc)) fail(`${solidName}.truncate(0.3) did not produce a valid mesh`);
        if (trunc) trunc.delete();

        // relax(int) + its C++-side clamp: relax(1e9) (INT32-valid) must clamp to
        // MAX_RELAX_ITERATIONS, not loop a billion times. Exercises the C++ clamp
        // only, not embind's double->int coercion near INT32_MAX.
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

  // Spline exports used by splines.html. Pin numeric behavior (not just finite
  // {x,y,z}) so a transposed-arg or wrong-target binding fails.
  const isVec = (v) => v && Number.isFinite(v.x) && Number.isFinite(v.y) && Number.isFinite(v.z);
  const approxVec = (v, x, y, z, eps = 1e-4) =>
    isVec(v) && Math.abs(v.x - x) <= eps && Math.abs(v.y - y) <= eps && Math.abs(v.z - z) <= eps;

  // p0..p3 control polygon (all unit vectors), shared by the cubic evaluators.
  const CTRL = [0, 1, 0, 1, 0, 0, 0, 0, 1, -1, 0, 0];
  for (const name of ['spline_cubic_fast', 'spline_cubic_slerp']) {
    const mid = Module[name](...CTRL, 0.5);
    if (!isVec(mid)) fail(`${name} returned non-finite ${JSON.stringify(mid)}`);
    // A cubic evaluator interpolates its endpoints: p0 at t=0, p3 at t=1 (both
    // unit, so normalized endpoints are exact). Catches arg transposition.
    const at0 = Module[name](...CTRL, 0);
    const at1 = Module[name](...CTRL, 1);
    if (!approxVec(at0, 0, 1, 0)) fail(`${name}(t=0) should be p0=(0,1,0), got ${JSON.stringify(at0)}`);
    if (!approxVec(at1, -1, 0, 0)) fail(`${name}(t=1) should be p3=(-1,0,0), got ${JSON.stringify(at1)}`);
  }
  // cubic_fast (normalized Bézier) and cubic_slerp (spherical de Casteljau)
  // differ away from the endpoints; identical midpoints mean both bindings
  // resolve to the same function (embind drift).
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
  // Pin cp1 = slerp(start, slerp(prev,end,.5), .5). Both slerps join orthonormal
  // vectors at t=0.5, so cp1 = normalize(start + normalize(prev+end)) =
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
