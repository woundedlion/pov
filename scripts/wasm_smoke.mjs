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

  // Run-wide: at least one pixel somewhere in the sweep must be non-zero. A
  // single effect legitimately renders black, but an all-zero sweep means the
  // draw path or the framebuffer view is dead — which the length check below
  // cannot see. Scanning stops once the flag is set, so this costs one partial
  // buffer walk for the whole run.
  let litPixel = false;

  // Run-wide counter for the per-frame telemetry the simulator reads: a single
  // effect reporting false is legitimate, an all-false sweep is not.
  let strobing = 0;

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
        // getBufferLength() describes the same span without touching the view.
        const len = engine.getBufferLength();
        if (len !== expected) {
          fail(`${name}: getBufferLength() ${len}, expected ${expected}`);
        }

        const strobe = engine.strobeColumns();
        if (typeof strobe !== 'boolean') {
          fail(`${name}: strobeColumns() returned ${typeof strobe} "${strobe}", expected a boolean`);
        } else if (strobe) {
          strobing++;
        }

        if (!litPixel) {
          for (let i = 0; i < px.length; i++) {
            if (px[i] !== 0) { litPixel = true; break; }
          }
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

    if (!litPixel) {
      fail('every effect at every resolution produced an all-zero framebuffer — ' +
        'the render path or getPixels() is not writing pixels');
    }
    if (strobing === 0) {
      fail('strobeColumns() reported false for every effect — the whole roster ' +
        'declares strobe = true');
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

    // ── Embind state seam: getParamGeneration / setAnimationsPaused ───────────
    {
      const effectNames = Object.keys(engine.getEffectSizes());
      if (effectNames.length === 0) {
        fail('state-seam: no effects to drive');
      } else {
        // The generation is the only token joining a definitions snapshot to a
        // later value read: it must hold across frames, reads and param writes,
        // hold across a rejected load, and advance by exactly one per accepted
        // load (effect_loads in wasm.cpp).
        const g0 = engine.getParamGeneration();
        if (!Number.isInteger(g0) || g0 < 1) {
          fail(`state-seam: getParamGeneration() = ${g0}, expected a positive integer ` +
            `after the sweep's effect loads`);
        }
        engine.drawFrame();
        engine.getParameterDefinitions();
        engine.getParamValues();
        engine.setParameter('definitely_not_a_param', 1);
        if (engine.getParamGeneration() !== g0) {
          fail(`state-seam: getParamGeneration() moved ${g0} -> ` +
            `${engine.getParamGeneration()} without an effect load`);
        }
        if (engine.setEffect('definitely_not_an_effect')) {
          fail('state-seam: setEffect(unknown) returned true');
        }
        if (engine.getParamGeneration() !== g0) {
          fail(`state-seam: a rejected setEffect bumped the generation to ` +
            `${engine.getParamGeneration()} (want ${g0})`);
        }
        for (let load = 1; load <= 2; load++) {
          if (!engine.setEffect(effectNames[0])) {
            fail(`state-seam: setEffect("${effectNames[0]}") failed`);
            break;
          }
          if (engine.getParamGeneration() !== g0 + load) {
            fail(`state-seam: generation ${engine.getParamGeneration()} after load ` +
              `${load}, expected ${g0 + load}`);
          }
        }

        if (RESOLUTIONS.length > 1) {
          const [w, h] = RESOLUTIONS[1];
          const beforeResolution = engine.getParamGeneration();
          if (!engine.setResolution(w, h)) {
            fail(`state-seam: setResolution(${w}, ${h}) failed`);
          } else {
            if (engine.getParamGeneration() !== beforeResolution + 1) {
              fail(`state-seam: resolution teardown left generation at ` +
                `${engine.getParamGeneration()}, expected ${beforeResolution + 1}`);
            }
            if (engine.getParameterDefinitions().length !== 0 ||
                engine.getParamValues().length !== 0) {
              fail('state-seam: resolution teardown left a non-empty parameter stream');
            }
            const noEffectGeneration = engine.getParamGeneration();
            engine.setResolution(w, h);
            if (engine.getParamGeneration() !== noEffectGeneration) {
              fail('state-seam: no-op setResolution changed the generation');
            }
            if (!engine.setEffect(effectNames[0])) {
              fail(`state-seam: setEffect("${effectNames[0]}") after resolution failed`);
            } else if (engine.getParamGeneration() !== noEffectGeneration + 1) {
              fail(`state-seam: post-resolution load generation ` +
                `${engine.getParamGeneration()}, expected ${noEffectGeneration + 1}`);
            }
          }
        }

        // setAnimationsPaused: pick the first effect whose animated params move
        // over a frame window, then assert the pause holds them still and the
        // resume lets them move again. animated=false params are engine-written
        // telemetry the pause gate does not cover.
        const animatedValues = () => engine.getParameterDefinitions()
          .filter((d) => d.animated).map((d) => d.value);
        const PAUSE_FRAMES = 20;
        const anyMoved = (before, after) => after.some((v, i) => v !== before[i]);
        let pauseTested = false;
        for (const name of effectNames) {
          if (!engine.setEffect(name)) continue;
          engine.drawFrame();
          const before = animatedValues();
          if (before.length === 0) continue;
          for (let f = 0; f < PAUSE_FRAMES; f++) engine.drawFrame();
          if (!anyMoved(before, animatedValues())) continue;

          engine.setAnimationsPaused(true);
          const held = animatedValues();
          for (let f = 0; f < PAUSE_FRAMES; f++) engine.drawFrame();
          const stillHeld = animatedValues();
          const moved = stillHeld.filter((v, i) => v !== held[i]).length;
          if (moved > 0) {
            fail(`state-seam: setAnimationsPaused(true) on ${name} left ${moved} ` +
              `animated param(s) moving`);
          }
          engine.setAnimationsPaused(false);
          for (let f = 0; f < PAUSE_FRAMES; f++) engine.drawFrame();
          if (!anyMoved(stillHeld, animatedValues())) {
            fail(`state-seam: setAnimationsPaused(false) on ${name} did not resume animation`);
          }
          console.log(`  state-seam: paramGeneration + setAnimationsPaused freeze/resume on ${name} OK`);
          pauseTested = true;
          break;
        }
        if (!pauseTested) {
          fail('state-seam: no effect had a moving animated param to exercise ' +
            'setAnimationsPaused');
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

      // ── getRecipe: recipe payload + reconstruction parity ──────────────────
      // For every entry with a recipe, replaying the chain from the seed through
      // the MeshOps op bindings must reproduce fromSolidName(entry) exactly —
      // same V/F/I counts and an identical vertex buffer (spec gate 8).
      if (typeof MeshOps.getRecipe !== 'function') {
        fail('MeshOps.getRecipe binding is missing');
      } else {
        // Unknown and recipe-less names return null, never abort.
        if (MeshOps.getRecipe('definitely_not_a_solid') !== null) {
          fail('getRecipe(unknown) should return null');
        }
        // Every islamic entry now carries a recipe mirror, so the recipe-less
        // cases live in the simple registry only.
        for (const bare of ['cube', 'octahedron']) {
          if (MeshOps.getRecipe(bare) !== null) {
            fail(`getRecipe("${bare}") should return null (no recipe)`);
          }
        }

        // Applies one recipe step through the wrapper's op bindings; params are
        // engine-native (radians for hankin, raw t, relax iteration count), so
        // they pass through verbatim.
        const applyStep = (wrapper, step) => {
          switch (step.op) {
            case 'kis': case 'ambo': case 'gyro': case 'dual':
            case 'meta': case 'needle': case 'zip':
              return wrapper[step.op]();
            case 'truncate': case 'bevel': case 'chamfer': case 'expand':
              return wrapper[step.op](step.param);
            case 'hankin': return wrapper.hankin(step.param);
            case 'relax': return wrapper.relax(step.param);
            case 'snub': return wrapper.snub(step.param, step.twist);
            default: return null;
          }
        };

        const RECIPE_ENTRIES = ['dodecahedron_hk62_ambo_hk62', 'octahedron_hk17_ambo_hk73'];
        for (const entryName of RECIPE_ENTRIES) {
          const recipe = MeshOps.getRecipe(entryName);
          if (!recipe) { fail(`getRecipe("${entryName}") returned null`); continue; }

          // Payload shape: seed string resolving in the registry, ops array of
          // {op: string, param: number, twist: number}.
          if (typeof recipe.seed !== 'string' || recipe.seed.length === 0) {
            fail(`${entryName}: recipe.seed "${recipe.seed}" is not a non-empty string`);
            continue;
          }
          if (!registry.some((e) => e.name === recipe.seed)) {
            fail(`${entryName}: recipe.seed "${recipe.seed}" not found in getRegistry()`);
            continue;
          }
          if (!Array.isArray(recipe.ops) || recipe.ops.length === 0) {
            fail(`${entryName}: recipe.ops is not a non-empty array`);
            continue;
          }
          let shapeOk = true;
          for (let i = 0; i < recipe.ops.length; i++) {
            const s = recipe.ops[i];
            if (typeof s.op !== 'string' || !Number.isFinite(s.param) || !Number.isFinite(s.twist)) {
              fail(`${entryName}: ops[${i}] malformed: ${JSON.stringify(s)}`);
              shapeOk = false;
            }
          }
          if (!shapeOk) continue;

          // Reconstruction parity against the whole-generator truth.
          let cur = MeshOps.fromSolidName(recipe.seed);
          if (!cur) { fail(`${entryName}: fromSolidName("${recipe.seed}") returned null`); continue; }
          let broken = false;
          for (let i = 0; i < recipe.ops.length; i++) {
            const next = applyStep(cur, recipe.ops[i]);
            cur.delete();
            cur = next;
            if (!cur) {
              fail(`${entryName}: ops[${i}] ${JSON.stringify(recipe.ops[i])} returned null`);
              broken = true;
              break;
            }
          }
          if (broken) continue;

          const truth = MeshOps.fromSolidName(entryName);
          if (!truth) {
            fail(`${entryName}: fromSolidName returned null`);
            cur.delete();
            continue;
          }
          const rv = cur.getVertices(), tv = truth.getVertices();
          const rf = cur.getFaces(), tf = truth.getFaces();
          if (rv.length !== tv.length) {
            fail(`${entryName}: reconstructed vertex count ${rv.length / 3} != ${tv.length / 3}`);
          } else if (rf.counts.length !== tf.counts.length) {
            fail(`${entryName}: reconstructed face count ${rf.counts.length} != ${tf.counts.length}`);
          } else if (rf.indices.length !== tf.indices.length) {
            fail(`${entryName}: reconstructed index count ${rf.indices.length} != ${tf.indices.length}`);
          } else {
            let mismatch = -1;
            for (let i = 0; i < rv.length; i++) {
              if (rv[i] !== tv[i]) { mismatch = i; break; }
            }
            if (mismatch >= 0) {
              fail(`${entryName}: reconstructed vertex buffer diverges at float ${mismatch}: ` +
                `${rv[mismatch]} != ${tv[mismatch]}`);
            }
          }
          cur.delete();
          truth.delete();
        }
        // The reconstructions leave finalized meshes in the tooling arena; wipe
        // so later sections start from a clean slate.
        MeshOps.clearToolingMemory();
        console.log(`  getRecipe: null cases + payload + reconstruction parity (${RECIPE_ENTRIES.join(', ')}) OK`);
      }

      // ── Operator roster: one topology signature per bound operator ──────────
      // MESHOP_BIND guarantees a .function() exists per MESHOP_LIST /
      // MESHOP_IRREGULAR_LIST name, not that it reaches the matching method —
      // two names on one method compile clean. Each row pins what the operator
      // is defined to produce from a cube (V=8, E=12, F=6): vertex/face/index
      // counts plus the face-degree and vertex-incidence histograms, per the
      // element formulas in core/mesh/conway.h (dual F/V swap, kis V+F and 2E
      // faces, ambo E and V+F, truncate 2E and V+F, expand 2E and V+E+F,
      // chamfer V+2E and F+E, snub 2E and V+2E+F) and the compositions
      // gyro = d·snub, meta = k·d·a, needle = k·d, zip = d·k, bevel = t·a.
      //
      // kis/needle and truncate/zip share all three counts on every seed, so
      // the histograms — not the counts — separate those pairs. The seed choice
      // is load-bearing too: on a self-dual seed truncate and zip coincide.
      const OP_SEED = 'cube';
      // "<degree>x<count>", ascending by degree.
      const histogram = (degrees) => {
        const h = new Map();
        for (const n of degrees) h.set(n, (h.get(n) ?? 0) + 1);
        return [...h.keys()].sort((a, b) => a - b).map((k) => `${k}x${h.get(k)}`).join(' ');
      };
      const signature = (w) => {
        const V = w.getVertices().length / 3;
        const f = w.getFaces();
        const incidence = new Array(V).fill(0);
        for (const i of f.indices) incidence[i]++;
        return `V=${V} F=${f.counts.length} I=${f.indices.length} ` +
          `faces[${histogram(f.counts)}] verts[${histogram(incidence)}]`;
      };

      // Fraction args stay at 0.25: truncate(0.5) merges each edge's two cut
      // vertices and collapses to ambo, a different topology.
      const FRACTION = 0.25;
      const OP_SIGNATURES = [
        ['dual', (s) => s.dual(), 'V=6 F=8 I=24 faces[3x8] verts[4x6]'],
        ['kis', (s) => s.kis(), 'V=14 F=24 I=72 faces[3x24] verts[4x6 6x8]'],
        ['ambo', (s) => s.ambo(), 'V=12 F=14 I=48 faces[3x8 4x6] verts[4x12]'],
        ['gyro', (s) => s.gyro(), 'V=38 F=24 I=120 faces[5x24] verts[3x32 4x6]'],
        ['meta', (s) => s.meta(), 'V=26 F=48 I=144 faces[3x48] verts[4x12 6x8 8x6]'],
        ['needle', (s) => s.needle(), 'V=14 F=24 I=72 faces[3x24] verts[3x8 8x6]'],
        ['zip', (s) => s.zip(), 'V=24 F=14 I=72 faces[4x6 6x8] verts[3x24]'],
        ['truncate', (s) => s.truncate(FRACTION), 'V=24 F=14 I=72 faces[3x8 8x6] verts[3x24]'],
        ['bevel', (s) => s.bevel(FRACTION), 'V=48 F=26 I=144 faces[4x12 6x8 8x6] verts[3x48]'],
        ['chamfer', (s) => s.chamfer(FRACTION), 'V=32 F=18 I=96 faces[4x6 6x12] verts[3x32]'],
        ['expand', (s) => s.expand(FRACTION), 'V=24 F=26 I=96 faces[3x8 4x18] verts[4x24]'],
        ['snub', (s) => s.snub(FRACTION, 0), 'V=24 F=38 I=120 faces[3x32 4x6] verts[5x24]'],
      ];

      // The table only catches an aliased name while its rows stay distinct.
      const distinct = new Set(OP_SIGNATURES.map(([, , want]) => want));
      if (distinct.size !== OP_SIGNATURES.length) {
        fail(`operator roster: ${OP_SIGNATURES.length} operators share ${distinct.size} ` +
          `expected signatures — the table cannot see a name aliased onto another method`);
      }

      const opSeed = MeshOps.fromSolidName(OP_SEED);
      if (!opSeed) {
        fail(`operator roster: fromSolidName("${OP_SEED}") returned null`);
      } else {
        const seedSig = signature(opSeed);
        const wantSeed = 'V=8 F=6 I=24 faces[4x6] verts[3x8]';
        if (seedSig !== wantSeed) {
          fail(`operator roster: seed ${OP_SEED} is ${seedSig}, expected ${wantSeed} — ` +
            `the expected operator signatures are derived from that seed`);
        } else {
          for (const [opName, apply, want] of OP_SIGNATURES) {
            const out = apply(opSeed);
            if (!out) {
              fail(`operator roster: ${OP_SEED}.${opName}() returned null`);
              continue;
            }
            const got = signature(out);
            if (got !== want) {
              fail(`operator roster: ${opName}(${OP_SEED}) is ${got}, expected ${want}`);
            }
            out.delete();
          }
          console.log(`  operator roster: ${OP_SIGNATURES.length} operator topologies pinned on ${OP_SEED}`);
        }
        opSeed.delete();
      }
      MeshOps.clearToolingMemory();
    }
  }

  // ── Color / palette / geometry exports ──────────────────────────────────────
  // These free functions and PaletteOps.bakeLut back the JS tool ports but are
  // never touched by the engine loop above; pin numeric behavior so a
  // transposed-arg or wrong-target binding fails here instead of shipping green.
  console.log('\nColor / palette / geometry:');

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

  // named_procedural_palettes: the browser tool's mirror of the palettes.h
  // X-macro. Shape and name uniqueness for every entry, plus two entries pinned
  // against the header literals; MAUVE_FADE's per-channel-asymmetric rows catch a
  // channel or an a/b/c/d transposition that a symmetric palette would hide.
  {
    const pals = Module.named_procedural_palettes();
    if (!Array.isArray(pals) || pals.length === 0) {
      fail(`named_procedural_palettes() returned ${JSON.stringify(pals)}, expected a non-empty array`);
    } else {
      const seen = new Set();
      for (const p of pals) {
        if (typeof p.name !== 'string' || p.name.length === 0) {
          fail(`named_procedural_palettes(): entry with no name: ${JSON.stringify(p)}`);
          continue;
        }
        if (seen.has(p.name)) fail(`named_procedural_palettes() repeats "${p.name}"`);
        seen.add(p.name);
        for (const k of ['a', 'b', 'c', 'd']) {
          const v = p[k];
          if (!Array.isArray(v) || v.length !== 3 || !v.every(Number.isFinite)) {
            fail(`${p.name}: coefficient ${k} = ${JSON.stringify(v)}, expected 3 finite floats`);
          }
        }
      }
      const PINNED = {
        DARK_RAINBOW: {
          a: [0.367, 0.367, 0.367], b: [0.5, 0.5, 0.5],
          c: [1, 1, 1], d: [0, 0.33, 0.67],
        },
        MAUVE_FADE: {
          a: [0.583, 0, 0.583], b: [1, 0, 1],
          c: [0.191, 0.348, 0.191], d: [0.175, 0.045, 0.15],
        },
      };
      for (const [name, want] of Object.entries(PINNED)) {
        const p = pals.find((e) => e.name === name);
        if (!p) {
          fail(`named_procedural_palettes() omits ${name}`);
          continue;
        }
        for (const k of ['a', 'b', 'c', 'd']) {
          const got = p[k];
          if (!Array.isArray(got) || got.some((v, i) => !approx(v, want[k][i], 1e-6))) {
            fail(`${name}.${k} = ${JSON.stringify(got)}, expected ${JSON.stringify(want[k])}`);
          }
        }
      }
      console.log(`  named_procedural_palettes: ${pals.length} entries, ` +
        `${Object.keys(PINNED).join(' + ')} coefficients pinned`);
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

  // Mobius: identity fixes every point; f(z) = 1/z is a 180° turn about x.
  {
    const v = { x: 0.48, y: 0.6, z: 0.64 };
    const id = Module.mobius_transform(v.x, v.y, v.z, 1, 0, 0, 0, 0, 0, 1, 0);
    if (!approxVec(id, v.x, v.y, v.z, 1e-3)) {
      fail(`mobius_transform identity = ${JSON.stringify(id)}, expected ${JSON.stringify(v)}`);
    }
    // a=0, b=1, c=1, d=0 in the eight-float order; a transposed binding breaks this.
    const inv = Module.mobius_transform(v.x, v.y, v.z, 0, 0, 1, 0, 1, 0, 0, 0);
    if (!approxVec(inv, v.x, -v.y, -v.z, 1e-3)) {
      fail(`mobius_transform 1/z = ${JSON.stringify(inv)}, expected (${v.x}, ${-v.y}, ${-v.z})`);
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
  console.log('  color/palette/geometry: transfer, interp, OKLab, HSV, procedural, lissajous, mobius, bakeLut OK');

  if (failures > 0) {
    console.error(`\nwasm_smoke: ${failures} failure(s)`);
    process.exitCode = 1;
    return;
  }
  console.log('\nwasm_smoke: OK');
}

await main();
