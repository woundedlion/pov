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

    // ── Embind write seam: setResolution / setClip / setParameter ─────────────
    // The per-effect loop above only READS the param streams; drive the write
    // methods end-to-end through embind so a binding-signature drift on the
    // setters fails here instead of shipping unseen.
    {
      const [w, h] = RESOLUTIONS[0];
      if (!engine.setResolution(w, h)) {
        fail(`write-seam: setResolution(${w}, ${h}) rejected a supported size`);
      }
      // Rejection path: an unsupported size returns false and keeps the prior
      // valid state (host predicate tests cover the predicate, not this seam).
      if (engine.setResolution(1, 1)) {
        fail('write-seam: setResolution(1, 1) accepted an unsupported size');
      }

      const effectNames = Object.keys(engine.getEffectSizes());
      if (effectNames.length === 0) {
        fail(`write-seam: no effects at ${w}x${h}`);
      } else {
        // setClip through embind: an in-range full-canvas band is accepted; a
        // negative, over-extent, or inverted band is rejected (false), never
        // trapped. The range check precedes the needs_full_frame branch, so both
        // outcomes are deterministic regardless of the effect.
        if (!engine.setEffect(effectNames[0])) {
          fail(`write-seam: setEffect("${effectNames[0]}") failed`);
        } else {
          if (!engine.setClip(0, w, 0, h)) fail('write-seam: setClip in-range full canvas rejected');
          if (engine.setClip(-1, w, 0, h)) fail('write-seam: setClip accepted a negative bound');
          if (engine.setClip(0, w + 1, 0, h)) fail('write-seam: setClip accepted x1 past the canvas width');
          if (engine.setClip(0, w, 0, h + 1)) fail('write-seam: setClip accepted y1 past the canvas height');
          if (engine.setClip(w, 0, 0, h)) fail('write-seam: setClip accepted an inverted (x0 > x1) band');
        }

        // setParameter unknown-name rejection.
        if (engine.setParameter('definitely_not_a_param', 1)) {
          fail('write-seam: setParameter(unknown name) returned true');
        }

        // setParameter clamp-readback: find a non-readonly float param with a
        // finite range, write past each bound, and read the effective value back
        // through getParameterDefinitions() (no drawFrame between, so animation
        // cannot move it). setParameter returns true even when it clamps, so the
        // readback — not the flag — is what proves the clamp took.
        let clampTested = false;
        const near = (a, b) => Number.isFinite(a) && Math.abs(a - b) <= 1e-3 * (1 + Math.abs(b));
        for (const name of effectNames) {
          if (!engine.setEffect(name)) continue;
          const defs = engine.getParameterDefinitions();
          if (!Array.isArray(defs)) continue;
          const t = defs.find((d) => typeof d.value === 'number' && !d.readonly &&
            Number.isFinite(d.min) && Number.isFinite(d.max) && d.max > d.min);
          if (!t) continue;
          const readBack = () => {
            const d = engine.getParameterDefinitions().find((e) => e.name === t.name);
            return d ? d.value : undefined;
          };
          const span = t.max - t.min;
          if (!engine.setParameter(t.name, t.max + span + 1)) {
            fail(`write-seam: setParameter("${t.name}", above max) was rejected`);
          } else if (!near(readBack(), t.max)) {
            fail(`write-seam: setParameter above max not clamped: read ${readBack()}, want ${t.max}`);
          }
          if (!engine.setParameter(t.name, t.min - span - 1)) {
            fail(`write-seam: setParameter("${t.name}", below min) was rejected`);
          } else if (!near(readBack(), t.min)) {
            fail(`write-seam: setParameter below min not clamped: read ${readBack()}, want ${t.min}`);
          }
          console.log(`  write-seam: setParameter clamp on ${name}."${t.name}" [${t.min}, ${t.max}] OK`);
          clampTested = true;
          break;
        }
        if (!clampTested) {
          fail('write-seam: found no float param to exercise the setParameter clamp');
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

  // ── MeshOps tooling bindings ────────────────────────────────────────────────
  // The engine loop drives only HolosphereEngine; exercise the MeshOpsWrapper
  // surface (used by solids.html) so its embind signatures can't drift unseen.
  console.log('\nMeshOps:');

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

  // ── Color / palette / lissajous exports ─────────────────────────────────────
  // These free functions and PaletteOps.bakeLut back the JS tool ports but are
  // never touched by the engine loop above; pin numeric behavior so a
  // transposed-arg or wrong-target binding fails here instead of shipping green.
  console.log('\nColor / palette / lissajous:');

  const approx = (a, b, eps = 1e-3) => Number.isFinite(a) && Math.abs(a - b) <= eps;
  const isVec = (v) => v && Number.isFinite(v.x) && Number.isFinite(v.y) && Number.isFinite(v.z);
  const approxVec = (v, x, y, z, eps = 1e-4) =>
    isVec(v) && Math.abs(v.x - x) <= eps && Math.abs(v.y - y) <= eps && Math.abs(v.z - z) <= eps;

  // sRGB transfer and its inverse: pinned endpoints plus a round-trip.
  {
    const s2l = Module.srgb_to_linear_float, l2s = Module.linear_to_srgb_float;
    if (!approx(s2l(0), 0) || !approx(s2l(1), 1)) fail(`srgb_to_linear_float endpoints off: ${s2l(0)}, ${s2l(1)}`);
    if (!approx(l2s(0), 0) || !approx(l2s(1), 1)) fail(`linear_to_srgb_float endpoints off: ${l2s(0)}, ${l2s(1)}`);
    const mid = s2l(0.5);
    if (!(mid > 0 && mid < 1)) fail(`srgb_to_linear_float(0.5) = ${mid}, expected interior (0,1)`);
    if (!approx(l2s(s2l(0.5)), 0.5, 2e-3)) fail(`sRGB transfer round-trip broken: ${l2s(s2l(0.5))}`);
  }

  // 16-bit linear interp LUT: 0 at black, saturating near 65535 at white, monotone.
  {
    const f0 = Module.srgb_to_linear_interp(0), f1 = Module.srgb_to_linear_interp(1),
      fm = Module.srgb_to_linear_interp(0.5);
    if (f0 !== 0) fail(`srgb_to_linear_interp(0) = ${f0}, expected 0`);
    if (!(f1 > 60000 && f1 <= 65535)) fail(`srgb_to_linear_interp(1) = ${f1}, expected near 65535`);
    if (!(fm > f0 && fm < f1)) fail(`srgb_to_linear_interp not monotone: ${f0} ${fm} ${f1}`);
  }

  // OKLab: white maps to L~1, a~0, b~0; an asymmetric color round-trips through both.
  {
    const white = Module.linear_rgb_to_oklab(1, 1, 1);
    if (!white || !approx(white.L, 1) || !approx(white.a, 0, 3e-3) || !approx(white.b, 0, 3e-3)) {
      fail(`linear_rgb_to_oklab(1,1,1) = ${JSON.stringify(white)}, expected ~{L:1,a:0,b:0}`);
    }
    // Round-trip an asymmetric linear color: catches r/g/b transposition either way.
    const lab = Module.linear_rgb_to_oklab(0.6, 0.3, 0.1);
    const rgb = lab && Module.oklab_to_linear_rgb(lab.L, lab.a, lab.b);
    if (!rgb || !approx(rgb.r, 0.6, 2e-3) || !approx(rgb.g, 0.3, 2e-3) || !approx(rgb.b, 0.1, 2e-3)) {
      fail(`OKLab round-trip broke: (0.6,0.3,0.1) -> ${JSON.stringify(rgb)}`);
    }
  }

  // HSV sextant path: s=0 is value-driven gray; a saturated hue-0 is red-dominant.
  {
    const gray = Module.hsv_to_rgb(0, 0, 200);
    if (!gray || gray.r !== gray.g || gray.g !== gray.b || !(gray.r > 150)) {
      fail(`hsv_to_rgb(0,0,200) = ${JSON.stringify(gray)}, expected value-driven gray`);
    }
    const black = Module.hsv_to_rgb(0, 0, 0);
    if (!black || black.r !== 0 || black.g !== 0 || black.b !== 0) {
      fail(`hsv_to_rgb(0,0,0) = ${JSON.stringify(black)}, expected (0,0,0)`);
    }
    const red = Module.hsv_to_rgb(0, 255, 255);
    if (!red || !(red.r > red.g && red.r > red.b)) {
      fail(`hsv_to_rgb(0,255,255) = ${JSON.stringify(red)}, expected red-dominant`);
    }
  }

  // Procedural cosine palette: an r-only cosine leaves g/b dark, pinning the channel target.
  {
    const c = Module.procedural_palette_linear(0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    if (!c || !(c.r > 0) || c.g !== 0 || c.b !== 0) {
      fail(`procedural_palette_linear r-only = ${JSON.stringify(c)}, expected r>0, g=b=0`);
    }
  }

  // Lissajous is unit-by-construction; t=0 is (0,1,0) for any m1,m2,a.
  {
    const p0 = Module.lissajous(3, 2, 0.5, 0);
    if (!approxVec(p0, 0, 1, 0)) fail(`lissajous(t=0) = ${JSON.stringify(p0)}, expected (0,1,0)`);
    // m1=1,m2=2,a=0,t=π/4 -> (cos π/4, 0, sin π/4); distinguishes m1/m2/a and target axes.
    const pq = Module.lissajous(1, 2, 0, Math.PI / 4);
    if (!approxVec(pq, Math.SQRT1_2, 0, Math.SQRT1_2)) {
      fail(`lissajous(1,2,0,π/4) = ${JSON.stringify(pq)}, expected (0.7071, 0, 0.7071)`);
    }
  }

  // PaletteOps.bakeLut: 256*3 sRGB bytes; a two-key gradient must vary end to end.
  {
    const po = new Module.PaletteOps();
    try {
      const lut = po.bakeLut(0, 0, 255, 255, 160, 255, 255, 160, 255, 255);
      if (!lut || lut.length !== 256 * 3) {
        fail(`bakeLut length ${lut && lut.length}, expected ${256 * 3}`);
      } else if (lut[0] === lut[765] && lut[1] === lut[766] && lut[2] === lut[767]) {
        fail(`bakeLut gradient is flat end-to-end: [${lut[0]},${lut[1]},${lut[2]}]`);
      }
    } finally {
      po.delete();
    }
  }
  console.log('  color/palette/lissajous: transfer, interp, OKLab, HSV, procedural, lissajous, bakeLut OK');

  if (failures > 0) {
    console.error(`\nwasm_smoke: ${failures} failure(s)`);
    process.exitCode = 1;
    return;
  }
  console.log('\nwasm_smoke: OK');
}

await main();
