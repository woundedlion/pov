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
import {
  darknessProblems,
  paramStreamProblems,
  stackCreepBudget,
} from './wasm_smoke_predicates.mjs';

const DEFAULT_JS = 'build/wasm-release/holosphere_wasm.js';
const jsArg = process.argv[2] || process.env.WASM_JS || DEFAULT_JS;
const jsPath = isAbsolute(jsArg) ? jsArg : join(process.cwd(), jsArg);

// The gate depth, spelled once here so CI and `just smoke` drive the same run:
// 120 frames reaches the late-lifecycle events (frame-48 ShapeShifter cut, arena
// compaction) a short run never hits. WASM_SMOKE_FRAMES shortens it for an
// ad-hoc local run.
const FRAMES_PER_EFFECT = Number(process.env.WASM_SMOKE_FRAMES ?? 120);

// The stack has no allocator trap and stack_high_water_mark() saturates at
// capacity, so `hwm > capacity` can never fire for it. Gate on an absolute byte
// budget instead (meaningful against any build's stack size); stackCreepBudget()
// mins it with the capacity fraction. Creep tripwire, not a bound: the HWM
// under-reports (see stack_high_water_mark() in wasm.cpp).
// The 2048 default is calibrated on the -O3 release build; -O0 debug frames run
// severalfold larger, so ci.yml overrides via WASM_SMOKE_STACK_CEILING.
const STACK_HWM_CEILING_BYTES = Number(process.env.WASM_SMOKE_STACK_CEILING ?? 2048);
// Shader's dynamic reference sweep peaks near 3.2 KB and retains half of
// the release stack as headroom under this effect-specific ratchet.
const SHADERBALL_STACK_HWM_CEILING_BYTES = 4096;

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

  let failures = 0;
  const fail = (msg) => { console.error(`  FAIL: ${msg}`); failures++; };

  // The tool pages' CSP grants 'wasm-unsafe-eval' but not 'unsafe-eval', which
  // holds only while the glue generates no code at runtime — so the probe spans
  // the whole run, not just module creation: it is armed before the glue's own
  // module evaluation and stays armed across every engine, MeshOps and free-
  // function call below (restored, and asserted, at the end of main). Both
  // routes to a code generator are counted: `Function` as constructor or
  // callee, and `eval` — replacing the global binding turns any glue `eval(s)`
  // into an indirect call through this proxy. A build that lost
  // -sDYNAMIC_EXECUTION=0 / -sEMBIND_AOT=1 fails here instead of throwing CSP
  // errors in the browser.
  const RealFunction = globalThis.Function;
  const realEval = globalThis.eval;
  let dynamicExecCalls = 0;
  let firstDynamicExec = null;
  const noteDynamicExec = (what) => {
    dynamicExecCalls++;
    firstDynamicExec ??= `${what}\n${new Error().stack}`;
  };
  globalThis.Function = new Proxy(RealFunction, {
    construct: (t, args, nt) => { noteDynamicExec('new Function()'); return Reflect.construct(t, args, nt); },
    apply: (t, self, args) => { noteDynamicExec('Function()'); return Reflect.apply(t, self, args); },
  });
  globalThis.eval = new Proxy(realEval, {
    apply: (t, self, args) => { noteDynamicExec('eval()'); return Reflect.apply(t, self, args); },
  });
  const disarmDynamicExecProbe = () => {
    globalThis.Function = RealFunction;
    globalThis.eval = realEval;
  };

  let Module;
  try {
    const { default: createHolosphereModule } = await import(pathToFileURL(jsPath));
    // Surface engine-side hs::log output and any abort() so a trap is visible in
    // the CI log rather than a bare non-zero exit.
    Module = await createHolosphereModule({
      print: (s) => console.log(`[wasm] ${s}`),
      printErr: (s) => console.error(`[wasm:err] ${s}`),
    });
  } catch (e) {
    disarmDynamicExecProbe();
    throw e;
  }

  // Per-effect-per-resolution darkness: an effect whose draw path regresses to
  // an all-zero framebuffer is invisible to a run-wide "something lit" flag,
  // and the length checks below cannot see it either. Keys are "Name@WxH".
  const darkEffects = new Set();
  const sweptEffects = new Set();

  // Run-wide counter for the per-frame telemetry the simulator reads: a single
  // effect reporting false is legitimate, an all-false sweep is not.
  let strobing = 0;

  // Enumerated from the module (generated from HS_RESOLUTIONS) so a new
  // resolution gets coverage without editing this file.
  const RESOLUTIONS = Module.HolosphereEngine.getSupportedResolutions();
  if (!RESOLUTIONS || RESOLUTIONS.length === 0) {
    console.error('wasm_smoke: getSupportedResolutions() returned no resolutions');
    process.exitCode = 1;
    return;
  }

  // setResolution/setEffect report their outcome as embind enums; pin the value
  // rosters so a dropped/renamed enumerator fails here. RESIZED and
  // ALREADY_ACTIVE are both successes (the requested size is active either
  // way), so the sweep accepts either.
  const RS = Module.ResolutionSetResult;
  const ES = Module.EffectSetResult;
  for (const outcome of ['RESIZED', 'ALREADY_ACTIVE', 'UNSUPPORTED']) {
    if (!RS || RS[outcome] === undefined) fail(`Module.ResolutionSetResult.${outcome} is not bound`);
  }
  for (const outcome of ['INSTALLED', 'UNKNOWN_EFFECT', 'UNSUPPORTED_RESOLUTION']) {
    if (!ES || ES[outcome] === undefined) fail(`Module.EffectSetResult.${outcome} is not bound`);
  }
  if (!RS || !ES) {
    console.error('wasm_smoke: result enums missing; every call below would misread');
    process.exitCode = 1;
    return;
  }
  const resolutionOk = (r) => r === RS.RESIZED || r === RS.ALREADY_ACTIVE;

  const engine = new Module.HolosphereEngine();
  try {
    for (const [w, h] of RESOLUTIONS) {
      if (!resolutionOk(engine.setResolution(w, h))) {
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
        sweptEffects.add(name);
        if (engine.setEffect(name) !== ES.INSTALLED) {
          fail(`setEffect("${name}") was rejected at ${w}x${h}`);
          continue;
        }
        if (name === 'Shader') {
          const presetCount = engine.getPresetCount();
          for (let preset = 0; preset < presetCount; preset++) {
            if (!engine.selectPreset(preset)) {
              fail(`Shader: selectPreset(${preset}) of ${presetCount} failed`);
              continue;
            }
            engine.drawFrame();
          }
          if (!engine.selectPreset(0)) {
            fail('Shader: could not restore preset 0 after the preset sweep');
          }
          engine.setAnimationsPaused(false);
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

        let lit = false;
        for (let i = 0; i < px.length; i++) {
          if (px[i] !== 0) { lit = true; break; }
        }
        if (!lit) darkEffects.add(`${name}@${w}x${h}`);

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
        const stackCeiling = name === 'Shader'
          ? Math.max(STACK_HWM_CEILING_BYTES, SHADERBALL_STACK_HWM_CEILING_BYTES)
          : STACK_HWM_CEILING_BYTES;
        const stackGate = stackCreepBudget(stack, stackCeiling);
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
        // The render path repaints the canary after every effect load, so the
        // mark above never reflects construction depth. init_high_water_mark is
        // latched at load time and survives the repaint; it is a running max, so
        // the first effect whose read crosses the budget is the one that did.
        if (stack && typeof stack.init_high_water_mark !== 'number') {
          fail(`${name}: getArenaMetrics().stack omits init_high_water_mark`);
        } else if (stack && stack.init_high_water_mark === 0) {
          fail(`${name}: init stack high-water mark is 0 (init canary tracking is broken)`);
        } else if (stack && stack.init_high_water_mark >= stackGate) {
          fail(`${name}: init-time stack high-water mark ${stack.init_high_water_mark} of ` +
            `${stack.capacity} bytes exceeds the ${stackGate}-byte creep budget — ` +
            `a stack-hungry effect constructor`);
        }

        // Exercise the embind param seam (getParameterDefinitions() +
        // getParamValues()) the GUI rides every frame: assert the two streams
        // stay zippable and well-formed.
        for (const problem of paramStreamProblems(engine.getParameterDefinitions(),
          engine.getParamValues())) {
          fail(`${name}: ${problem}`);
        }
      }
    }

    // RingShower expands its rings from zero radius and lights around frame 24,
    // so a short window is black by design; tests/test_effects.h carries the
    // same exemption, as a static_assert on the roster entry.
    for (const problem of darknessProblems({
      frames: FRAMES_PER_EFFECT,
      resolutions: RESOLUTIONS,
      sweptEffects,
      darkKeys: darkEffects,
    })) {
      fail(problem);
    }
    console.log(`\nall-black: ${darkEffects.size} of ${sweptEffects.size * RESOLUTIONS.length} ` +
      `effect/resolution passes rendered no lit pixel` +
      (darkEffects.size ? ` (${[...darkEffects].sort().join(', ')})` : ''));

    if (strobing === 0) {
      fail('strobeColumns() reported false for every effect — the whole roster ' +
        'declares strobe = false, or the binding is dead');
    }

    // ── Embind write seam: setResolution / setClip / setParameter ─────────────
    // The per-effect loop above only READS the param streams; drive the write
    // methods end-to-end through embind so a binding-signature drift on the
    // setters fails here instead of shipping unseen.
    {
      const [w, h] = RESOLUTIONS[0];
      if (!resolutionOk(engine.setResolution(w, h))) {
        fail(`write-seam: setResolution(${w}, ${h}) rejected a supported size`);
      }
      // A repeat of the active size is the distinguished no-op, not a resize.
      if (engine.setResolution(w, h) !== RS.ALREADY_ACTIVE) {
        fail('write-seam: repeated setResolution did not report ALREADY_ACTIVE');
      }
      // Rejection path: an unsupported size reports UNSUPPORTED and keeps the
      // prior valid state (host predicate tests cover the predicate, not this
      // seam).
      if (engine.setResolution(1, 1) !== RS.UNSUPPORTED) {
        fail('write-seam: setResolution(1, 1) did not report UNSUPPORTED');
      }

      const effectNames = Object.keys(engine.getEffectSizes());
      if (effectNames.length === 0) {
        fail(`write-seam: no effects at ${w}x${h}`);
      } else {
        // setClip through embind: an in-range full-canvas band succeeds as
        // APPLIED, or as FULL_FRAME_KEPT when the effect reports
        // needs_full_frame(); effectNames[0] is not pinned to either, so accept
        // both. A negative, over-extent, or inverted band reports
        // INVALID_BOUNDS, never traps. The range check precedes the
        // needs_full_frame branch, so the rejection is deterministic regardless
        // of the effect. NaN, fractional and past-i32 numbers are the coercion
        // cases an int parameter would have delivered as 0 (an empty band under
        // an APPLIED result). INVALID_BOUNDS must stay distinct from NO_EFFECT
        // — the pool faults on the former and must not fault on the latter — so
        // pin the roster too.
        const C = Module.ClipSetResult;
        for (const outcome of ['APPLIED', 'NO_EFFECT', 'INVALID_BOUNDS', 'FULL_FRAME_KEPT']) {
          if (!C || C[outcome] === undefined) fail(`Module.ClipSetResult.${outcome} is not bound`);
        }
        if (engine.setEffect(effectNames[0]) !== ES.INSTALLED) {
          fail(`write-seam: setEffect("${effectNames[0]}") failed`);
        } else {
          const full = engine.setClip(0, w, 0, h);
          if (full !== C.APPLIED && full !== C.FULL_FRAME_KEPT) {
            fail('write-seam: setClip in-range full canvas rejected');
          }
          for (const [what, bounds] of [
            ['a negative bound', [-1, w, 0, h]],
            ['x1 past the canvas width', [0, w + 1, 0, h]],
            ['y1 past the canvas height', [0, w, 0, h + 1]],
            ['an inverted (x0 > x1) band', [w, 0, 0, h]],
            ['a NaN bound', [NaN, w, 0, h]],
            ['a fractional bound', [0, w - 0.5, 0, h]],
            ['a bound past the i32 range', [0, 2 ** 32, 0, h]],
          ]) {
            if (engine.setClip(...bounds) !== C.INVALID_BOUNDS) {
              fail(`write-seam: setClip did not report INVALID_BOUNDS for ${what}`);
            }
          }
          // A resolution change tears the effect down, so the clip that follows
          // reports NO_EFFECT — a benign ordering state, not a bad-bounds fault.
          const [w2, h2] = RESOLUTIONS[RESOLUTIONS.length - 1];
          if (w2 !== w || h2 !== h) {
            if (engine.setResolution(w2, h2) !== RS.RESIZED) fail(`write-seam: setResolution(${w2}, ${h2}) did not report RESIZED`);
            if (engine.setClip(0, w2, 0, h2) !== C.NO_EFFECT) {
              fail('write-seam: setClip with no effect set did not report NO_EFFECT');
            }
            if (engine.setResolution(w, h) !== RS.RESIZED) fail(`write-seam: setResolution(${w}, ${h}) did not report RESIZED`);
            if (engine.setEffect(effectNames[0]) !== ES.INSTALLED) {
              fail(`write-seam: setEffect("${effectNames[0]}") failed after the resolution round-trip`);
            }
          }
        }

        // setParameter reports its outcome as the ParamSetResult embind enum;
        // pin the value roster so a dropped/renamed enumerator fails here.
        const R = Module.ParamSetResult;
        for (const reason of
             ['APPLIED', 'NO_EFFECT', 'UNKNOWN_PARAM', 'READONLY', 'NON_FINITE']) {
          if (!R || !R[reason]) {
            fail(`write-seam: Module.ParamSetResult.${reason} is not bound`);
          }
        }

        // setParameter unknown-name rejection.
        if (engine.setParameter('definitely_not_a_param', 1) !== R.UNKNOWN_PARAM) {
          fail('write-seam: setParameter(unknown name) did not report UNKNOWN_PARAM');
        }

        // setParameter clamp-readback: find a non-readonly float param with a
        // finite range, write past each bound, and read the effective value back
        // through getParameterDefinitions() (no drawFrame between, so animation
        // cannot move it). setParameter reports APPLIED even when it clamps, so
        // the readback — not the result — is what proves the clamp took.
        let clampTested = false;
        const near = (a, b) => Number.isFinite(a) && Math.abs(a - b) <= 1e-3 * (1 + Math.abs(b));
        for (const name of effectNames) {
          if (engine.setEffect(name) !== ES.INSTALLED) continue;
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
          if (engine.setParameter(t.name, t.max + span + 1) !== R.APPLIED) {
            fail(`write-seam: setParameter("${t.name}", above max) was rejected`);
          } else if (!near(readBack(), t.max)) {
            fail(`write-seam: setParameter above max not clamped: read ${readBack()}, want ${t.max}`);
          }
          if (engine.setParameter(t.name, t.min - span - 1) !== R.APPLIED) {
            fail(`write-seam: setParameter("${t.name}", below min) was rejected`);
          } else if (!near(readBack(), t.min)) {
            fail(`write-seam: setParameter below min not clamped: read ${readBack()}, want ${t.min}`);
          }
          // A non-finite write names its reason and leaves the value in place.
          if (engine.setParameter(t.name, NaN) !== R.NON_FINITE) {
            fail(`write-seam: setParameter("${t.name}", NaN) did not report NON_FINITE`);
          } else if (!near(readBack(), t.min)) {
            fail(`write-seam: setParameter(NaN) moved the value to ${readBack()}`);
          }
          console.log(`  write-seam: setParameter clamp on ${name}."${t.name}" [${t.min}, ${t.max}] OK`);
          clampTested = true;
          break;
        }
        if (!clampTested) {
          fail('write-seam: found no float param to exercise the setParameter clamp');
        }
        // A readonly param write names its reason; skip if the roster exposes
        // none (engine_contract_wasm.test.js requires at least one, so drift
        // toward zero readonly params is caught there).
        let readonlyTested = false;
        for (const name of effectNames) {
          if (engine.setEffect(name) !== ES.INSTALLED) continue;
          const t = engine.getParameterDefinitions().find((d) => d.readonly);
          if (!t) continue;
          if (engine.setParameter(t.name, t.value) !== R.READONLY) {
            fail(`write-seam: setParameter(readonly "${t.name}") did not report READONLY`);
          }
          readonlyTested = true;
          break;
        }
        if (!readonlyTested) {
          console.log('  write-seam: no readonly param in the roster; READONLY untested');
        }
      }
      // The clamp probe picks its target by range, not by animation, and an
      // APPLIED write to an animated param latches the engine-wide pause, which
      // survives setEffect. Clear it so the state seam below starts unpaused.
      engine.setAnimationsPaused(false);
    }

    // ── Shader authoring route ──────────────────────────────────────────────
    {
      const R = Module.ParamSetResult;
      const close = (a, b) => Number.isFinite(a) && Number.isFinite(b) &&
        Math.abs(a - b) <= 1e-3 * (1 + Math.abs(b));
      if (engine.setEffect('ShaderBall') !== ES.INSTALLED) {
        fail('shader-authoring: legacy setEffect("ShaderBall") alias failed');
      }
      if (engine.setEffect('Shader') !== ES.INSTALLED) {
        fail('shader-authoring: setEffect("Shader") failed');
      } else {
        if (engine.setParameter('Lens', 2) !== R.APPLIED) {
          fail('shader-authoring: Twist lens write failed');
        }
        engine.drawFrame();
        let snapshot = engine.getFullConfigSnapshot();
        if (!snapshot || snapshot.pendingFieldIds.length !== 0) {
          fail('shader-authoring: valid uncompiled Twist lens stayed pending');
        }

        engine.setEffect('Shader');
        engine.setParameter('Planar Warp 1', 4);
        engine.setParameter('Planar Warp 1 Scale', 1);
        engine.setParameter('Planar Warp 1 Strength', 1);
        engine.setParameter('Planar Warp 1', 5);
        engine.drawFrame();
        const strength = engine.getParameterDefinitions().find(
          (definition) => definition.name === 'Planar Warp 1 Strength');
        if (!strength || !close(strength.value, strength.max)) {
          fail(`shader-authoring: Curl Flow strength ${strength?.value} ` +
            `was not clamped to ${strength?.max}`);
        }
        snapshot = engine.getFullConfigSnapshot();
        if (!snapshot || snapshot.pendingFieldIds.length !== 0) {
          fail('shader-authoring: valid Curl Flow configuration stayed pending');
        }

        engine.setParameter('Function', 6);
        engine.drawFrame();
        snapshot = engine.getFullConfigSnapshot();
        if (!snapshot || snapshot.pendingFieldIds.length === 0) {
          fail('shader-authoring: incompatible sphere source bypassed admission');
        }
      }
      engine.setAnimationsPaused(false);
    }

    // ── ShaderChain authoring route ─────────────────────────────────────────
    // setShaderChain / getShaderChainCatalog: the chain interpreter's embind
    // seam. The catalog is a class function (no engine needed); setShaderChain
    // compiles synchronously — on APPLIED the param definitions are already
    // rebuilt and the generation bumped before it returns.
    {
      let catalog = null;
      try {
        catalog = JSON.parse(Module.HolosphereEngine.getShaderChainCatalog());
      } catch {
        fail('shader-chain: getShaderChainCatalog() is not valid JSON');
      }
      if (catalog) {
        if (catalog.catalog_version !== 2) {
          fail(`shader-chain: catalog_version ${catalog.catalog_version}, expected 2`);
        }
        for (const key of ['max_chain_ops', 'arena_bytes', 'max_params',
          'max_instance_id_length', 'per_op_overhead_bytes',
          'per_param_name_bytes']) {
          if (!Number.isInteger(catalog.budgets?.[key]) || catalog.budgets[key] <= 0) {
            fail(`shader-chain: catalog budgets.${key} = ${catalog.budgets?.[key]}, ` +
              `expected a positive integer`);
          }
        }
        if (!Array.isArray(catalog.operators) || catalog.operators.length < 4) {
          fail(`shader-chain: catalog lists ${catalog.operators?.length} operators, expected >= 4`);
        }
      }

      const DEFAULT_CHAIN = [
        { instance: 'camera', operator: 'sphere.rotate.v2' },
        { instance: 'project', operator: 'project.stereographic.v2' },
        { instance: 'sample', operator: 'sample.grid.v2' },
        { instance: 'colorize', operator: 'colorize.generated-palette.v2' },
      ];

      // The refusal for a non-chain effect names itself, never traps.
      if (engine.setEffect('Shader') !== ES.INSTALLED) {
        fail('shader-chain: setEffect("Shader") failed');
      } else if (engine.setShaderChain(DEFAULT_CHAIN).code !== 'NOT_CHAIN_EFFECT') {
        fail('shader-chain: setShaderChain on a non-chain effect did not report NOT_CHAIN_EFFECT');
      }

      if (engine.setEffect('ShaderChain') !== ES.INSTALLED) {
        fail('shader-chain: setEffect("ShaderChain") failed');
      } else {
        const applied = engine.setShaderChain(DEFAULT_CHAIN);
        if (applied.code !== 'APPLIED' || applied.entryIndex !== -1) {
          fail(`shader-chain: default chain refused: ${applied.code}@${applied.entryIndex}`);
        }
        engine.drawFrame();
        const px = engine.getPixels();
        let lit = false;
        for (let i = 0; i < px.length; i++) {
          if (px[i] !== 0) { lit = true; break; }
        }
        if (!lit) fail('shader-chain: default chain rendered an all-black frame');

        const names = new Set(engine.getParameterDefinitions().map((d) => d.name));
        for (const want of ['camera.wander', 'project.pole-fade',
          'sample.pattern-freq', 'sample.coverage-mode', 'colorize.palette-mode']) {
          if (!names.has(want)) fail(`shader-chain: definitions omit "${want}"`);
        }
        const coverage = engine.getParameterDefinitions().find(
          (d) => d.name === 'sample.coverage-mode');
        if (!coverage || !Array.isArray(coverage.options) ||
            coverage.options.length !== 4) {
          fail('shader-chain: sample.coverage-mode did not surface its 4 enum options');
        }

        // A re-chain bumps the generation and swaps the instance namespace.
        const g0 = engine.getParamGeneration();
        const rechain = engine.setShaderChain([
          { instance: 'camera2', operator: 'sphere.rotate.v2' },
          ...DEFAULT_CHAIN.slice(1),
        ]);
        if (rechain.code !== 'APPLIED') {
          fail(`shader-chain: re-chain refused: ${rechain.code}`);
        }
        if (engine.getParamGeneration() === g0) {
          fail('shader-chain: an APPLIED re-chain did not bump getParamGeneration()');
        }
        const renamed = new Set(engine.getParameterDefinitions().map((d) => d.name));
        if (!renamed.has('camera2.wander') || renamed.has('camera.wander')) {
          fail('shader-chain: re-chain did not swap the instance parameter namespace');
        }

        // A refusal names its code and entry and leaves the prior state whole.
        const g1 = engine.getParamGeneration();
        const refused = engine.setShaderChain([
          { instance: 'camera', operator: 'sphere.rotate.v999' },
          ...DEFAULT_CHAIN.slice(1),
        ]);
        if (refused.code !== 'UNKNOWN_OPERATOR' || refused.entryIndex !== 0) {
          fail(`shader-chain: unknown operator reported ${refused.code}@${refused.entryIndex}`);
        }
        if (engine.getParamGeneration() !== g1) {
          fail('shader-chain: a refused chain moved getParamGeneration()');
        }
        if (!new Set(engine.getParameterDefinitions().map((d) => d.name)).has('camera2.wander')) {
          fail('shader-chain: a refused chain disturbed the registered definitions');
        }
        engine.drawFrame();
        const after = engine.getPixels();
        let stillLit = false;
        for (let i = 0; i < after.length; i++) {
          if (after[i] !== 0) { stillLit = true; break; }
        }
        if (!stillLit) fail('shader-chain: render went black after a refused chain');

        // The boundary rejects malformed payloads without trapping.
        for (const [what, payload] of [
          ['a non-array payload', 42],
          ['a null entry', [null]],
          ['a non-string instance', [{ instance: 7, operator: 'sphere.rotate.v2' }]],
          ['a missing operator field', [{ instance: 'camera' }]],
        ]) {
          if (engine.setShaderChain(payload).code !== 'MALFORMED_PAYLOAD') {
            fail(`shader-chain: ${what} did not report MALFORMED_PAYLOAD`);
          }
        }
        console.log('  shader-chain: catalog + setShaderChain apply/rebind/refusal/malformed OK');
      }
      engine.setAnimationsPaused(false);
    }

    // Stable fixed-effect and preset identities are independent of C++ names
    // and generated array positions.
    {
      if (engine.setEffect('mobius-grid') !== ES.INSTALLED) {
        fail('stable-identity: setEffect("mobius-grid") failed');
      } else {
        const ids = Array.from(engine.getPresetIds());
        if (ids.join(',') !== 'mobius-grid,mobius-grid-2') {
          fail(`stable-identity: unexpected Mobius preset IDs ${ids}`);
        }
        if (!engine.selectPresetById('mobius-grid-2') ||
            engine.getPresetIndex() !== 1) {
          fail('stable-identity: mobius-grid-2 selection failed');
        }
        if (engine.selectPresetById('missing-preset')) {
          fail('stable-identity: unknown preset was accepted');
        }
      }
      engine.setAnimationsPaused(false);
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
        // load.
        // The authoring probe above may leave a Shader schema refresh for
        // its next frame. Settle it before measuring generation immutability.
        engine.drawFrame();
        engine.getParameterDefinitions();
        engine.getParamValues();
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
        if (engine.setEffect('definitely_not_an_effect') !== ES.UNKNOWN_EFFECT) {
          fail('state-seam: setEffect(unknown) did not report UNKNOWN_EFFECT');
        }
        if (engine.getParamGeneration() !== g0) {
          fail(`state-seam: a rejected setEffect bumped the generation to ` +
            `${engine.getParamGeneration()} (want ${g0})`);
        }
        for (let load = 1; load <= 2; load++) {
          if (engine.setEffect(effectNames[0]) !== ES.INSTALLED) {
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
          if (engine.setResolution(w, h) !== RS.RESIZED) {
            fail(`state-seam: setResolution(${w}, ${h}) did not report RESIZED`);
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
            if (engine.setResolution(w, h) !== RS.ALREADY_ACTIVE) {
              fail('state-seam: no-op setResolution did not report ALREADY_ACTIVE');
            }
            if (engine.getParamGeneration() !== noEffectGeneration) {
              fail('state-seam: no-op setResolution changed the generation');
            }
            if (engine.setEffect(effectNames[0]) !== ES.INSTALLED) {
              fail(`state-seam: setEffect("${effectNames[0]}") after resolution failed`);
            } else if (engine.getParamGeneration() !== noEffectGeneration + 1) {
              fail(`state-seam: post-resolution load generation ` +
                `${engine.getParamGeneration()}, expected ${noEffectGeneration + 1}`);
            }
          }
        }

        // Pick the first effect whose animated params move over a frame window.
        // A manual animated write must engage pause, survive an effect reload,
        // and resume only after setAnimationsPaused(false).
        const animatedValues = () => engine.getParameterDefinitions()
          .filter((d) => d.animated).map((d) => d.value);
        const PAUSE_FRAMES = 20;
        const anyMoved = (before, after) => after.some((v, i) => v !== before[i]);
        let pauseTested = false;
        for (const name of effectNames) {
          if (engine.setEffect(name) !== ES.INSTALLED) continue;
          engine.drawFrame();
          const before = animatedValues();
          if (before.length === 0) continue;
          for (let f = 0; f < PAUSE_FRAMES; f++) engine.drawFrame();
          if (!anyMoved(before, animatedValues())) continue;

          const manual = engine.getParameterDefinitions().find((d) => d.animated);
          if (engine.setParameter(manual.name, manual.value) !==
              Module.ParamSetResult.APPLIED) {
            fail(`state-seam: manual animated write on ${name}.${manual.name} failed`);
          }
          if (engine.getAnimationsPaused() !== true) {
            fail(`state-seam: an animated write on ${name}.${manual.name} left ` +
              'getAnimationsPaused() false — the implicit pause is unreadable');
          }
          const held = animatedValues();
          for (let f = 0; f < PAUSE_FRAMES; f++) engine.drawFrame();
          const stillHeld = animatedValues();
          const moved = stillHeld.filter((v, i) => v !== held[i]).length;
          if (moved > 0) {
            fail(`state-seam: setAnimationsPaused(true) on ${name} left ${moved} ` +
              `animated param(s) moving`);
          }

          if (engine.setEffect(name) !== ES.INSTALLED) {
            fail(`state-seam: setEffect("${name}") while paused failed`);
            break;
          }
          const reloadedHeld = animatedValues();
          for (let f = 0; f < PAUSE_FRAMES; f++) engine.drawFrame();
          if (anyMoved(reloadedHeld, animatedValues())) {
            fail(`state-seam: setEffect("${name}") lost the paused state`);
          }
          if (engine.getAnimationsPaused() !== true) {
            fail(`state-seam: setEffect("${name}") left getAnimationsPaused() false ` +
              'while the effect was still frozen');
          }

          engine.setAnimationsPaused(false);
          for (let f = 0; f < PAUSE_FRAMES; f++) engine.drawFrame();
          if (!anyMoved(reloadedHeld, animatedValues())) {
            fail(`state-seam: setAnimationsPaused(false) on ${name} did not resume animation`);
          }
          if (engine.getAnimationsPaused() !== false) {
            fail(`state-seam: setAnimationsPaused(false) on ${name} left ` +
              'getAnimationsPaused() true');
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

    // ── Pole LOD seam: setPoleLod / getPoleLod ────────────────────────────────
    // The sweep above renders at the shipped default (0, decimation off), so
    // nothing else drives this pair. Write a non-zero aggressiveness, read it
    // back, and render through the decimated scan.
    {
      const effectNames = Object.keys(engine.getEffectSizes());
      if (effectNames.length === 0) {
        fail('pole-lod: no effects to drive');
      } else if (engine.setEffect(effectNames[0]) !== ES.INSTALLED) {
        fail(`pole-lod: setEffect("${effectNames[0]}") failed`);
      } else {
        const original = engine.getPoleLod();
        if (!Number.isFinite(original) || original < 0) {
          fail(`pole-lod: getPoleLod() = ${original}, expected a finite non-negative default`);
        }
        try {
          // 1.5 is exact in float32, so the readback needs no epsilon.
          engine.setPoleLod(1.5);
          if (engine.getPoleLod() !== 1.5) {
            fail(`pole-lod: getPoleLod() = ${engine.getPoleLod()} after setPoleLod(1.5)`);
          }
          engine.drawFrame();
          if (engine.getPixels().length === 0) {
            fail('pole-lod: getPixels() is empty after a decimated frame');
          }
          // Out-of-domain inputs clamp, never trap (setPoleLod in wasm.cpp).
          engine.setPoleLod(-1);
          if (engine.getPoleLod() !== 0) {
            fail(`pole-lod: setPoleLod(-1) left ${engine.getPoleLod()}, expected the 0 clamp`);
          }
          engine.setPoleLod(1e9);
          if (engine.getPoleLod() !== 8) {
            fail(`pole-lod: setPoleLod(1e9) left ${engine.getPoleLod()}, expected the 8 ceiling`);
          }
        } finally {
          // Later steps must see the value the sweep ran at.
          engine.setPoleLod(original);
        }
        console.log(`  pole-lod: setPoleLod/getPoleLod readback + decimated frame on ${effectNames[0]} OK`);
      }
    }

    // ── Memory-growth seam: a held view detaches, the re-fetch is live ────────
    // getPixels() aliases WASM linear memory and ALLOW_MEMORY_GROWTH detaches
    // that view on any heap growth — the bridge's most documented invariant, and
    // otherwise never driven end to end: the MeshOps section below runs after
    // engine.delete(), so no pixel view is ever outstanding when the growth
    // happens. MeshOps' 16 MB tooling block is allocated lazily on first use and
    // cannot fit the heap this module starts with, so the first call here IS the
    // growth. Must stay ahead of every other MeshOps call in this file.
    {
      const effectNames = Object.keys(engine.getEffectSizes());
      const registry = Module.MeshOps && Module.MeshOps.getRegistry();
      const growSolid = registry && registry.length ? registry[0].name : null;
      if (effectNames.length === 0) {
        fail('growth-seam: no effects to render a frame with');
      } else if (!growSolid) {
        fail('growth-seam: MeshOps reported no solid to allocate the tooling block with');
      } else if (engine.setEffect(effectNames[0]) !== ES.INSTALLED) {
        fail(`growth-seam: setEffect("${effectNames[0]}") failed`);
      } else {
        engine.drawFrame();
        const held = engine.getPixels();
        const before = Array.from(held);
        if (held.byteLength === 0) {
          fail('growth-seam: the freshly fetched pixel view is already detached');
        }
        const grower = Module.MeshOps.fromSolidName(growSolid);
        if (!grower) fail(`growth-seam: fromSolidName("${growSolid}") returned null`);
        else grower.delete();
        // Leave the arenas as the MeshOps section below expects to find them.
        Module.MeshOps.clearToolingMemory();

        if (held.byteLength !== 0) {
          fail('growth-seam: the 16 MB tooling allocation left the held pixel view ' +
            'attached — the detachment contract is no longer exercised, so a ' +
            'regression in it would ship unseen');
        }
        const refetched = engine.getPixels();
        if (refetched.byteLength === 0) {
          fail('growth-seam: the re-fetched pixel view is detached');
        } else if (refetched.length !== before.length) {
          fail(`growth-seam: re-fetched view length ${refetched.length}, ` +
            `expected ${before.length}`);
        } else {
          let mismatch = -1;
          for (let i = 0; i < before.length; i++) {
            if (refetched[i] !== before[i]) { mismatch = i; break; }
          }
          if (mismatch >= 0) {
            fail(`growth-seam: re-fetched pixels diverge at index ${mismatch}: ` +
              `${refetched[mismatch]} != ${before[mismatch]}`);
          } else {
            console.log(`  growth-seam: held view detached by the tooling allocation; ` +
              `re-fetch is live and byte-identical (${before.length} channels)`);
          }
        }
      }
    }

    // Log worst-case stack usage as a margin against STACK_SIZE.
    const stack = engine.getArenaMetrics().stack;
    if (!stack) fail('getArenaMetrics() omits the stack region');
    else console.log(`\nstack: ${stack.high_water_mark}/${stack.capacity} bytes peak ` +
      `(effect init peak ${stack.init_high_water_mark})`);
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
    // A rejected mesh-producing call returns null and names the reason through
    // the MeshOpResult embind enum; pin the value roster so a dropped/renamed
    // enumerator fails here. The reasons drive opposite caller actions (shrink
    // the chain vs. clearToolingMemory()), so an undifferentiated null is not
    // enough for the solids tool.
    const MR = Module.MeshOpResult;
    for (const reason of
         ['OK', 'UNKNOWN_NAME', 'CONNECTIVITY_OVERFLOW', 'FACE_DEGREE_OVERFLOW',
          'ARENA_EXHAUSTED', 'NON_FINITE_ARG', 'ANGLE_OUT_OF_DOMAIN',
          'STALE_WRAPPER']) {
      if (!MR || MR[reason] === undefined) {
        fail(`Module.MeshOpResult.${reason} is not bound`);
      }
    }
    if (typeof MeshOps.getLastResult !== 'function') {
      fail('MeshOps.getLastResult binding is missing');
    }

    // Use a real registry name rather than hardcoding one (anti-drift).
    const registry = MeshOps.getRegistry();
    const solidName = registry && registry.length ? registry[0].name : null;
    if (!solidName) {
      fail('MeshOps.getRegistry() returned no solids');
    } else {
      // Unknown names must be rejected (null), not abort the module.
      const bogus = MeshOps.fromSolidName('definitely_not_a_solid');
      if (bogus) { fail('fromSolidName(unknown) should return null'); bogus.delete(); }
      if (MeshOps.getLastResult() !== MR.UNKNOWN_NAME) {
        fail('fromSolidName(unknown) did not report UNKNOWN_NAME');
      }

      const solid = MeshOps.fromSolidName(solidName);
      if (MeshOps.getLastResult() !== MR.OK) {
        fail(`fromSolidName("${solidName}") did not report OK`);
      }
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
        if (MeshOps.getLastResult() !== MR.NON_FINITE_ARG) {
          fail(`${solidName}.hankin(NaN) did not report NON_FINITE_ARG`);
        }
        // A finite angle outside [0, pi/2] is a different rejection: it is
        // rejected rather than clamped, and must not read back as NON_FINITE_ARG.
        const hankinWide = solid.hankin(Math.PI);
        if (hankinWide) { fail(`${solidName}.hankin(π) should return null`); hankinWide.delete(); }
        if (MeshOps.getLastResult() !== MR.ANGLE_OUT_OF_DOMAIN) {
          fail(`${solidName}.hankin(π) did not report ANGLE_OUT_OF_DOMAIN`);
        }

        solid.delete();
      }

      // Tooling arena high-water marks must stay within capacity after the ops.
      // The two scratch arenas are what TOOLING_BYTES_PER_MESH_ELEMENT is sized
      // against, so a missing or unbound region leaves the loop below blind to
      // the regions an operator overrun would kill the module through.
      const tm = MeshOps.getArenaMetrics();
      for (const region of ['tooling_arena', 'tooling_scratch_a', 'tooling_scratch_b']) {
        if (!tm[region]) {
          fail(`MeshOps.getArenaMetrics() omits ${region}`);
        } else if (!(tm[region].capacity > 0)) {
          fail(`MeshOps ${region} capacity is ${tm[region].capacity} after the ops ran`);
        }
      }
      for (const region of Object.keys(tm)) {
        const { high_water_mark: hwm, capacity } = tm[region];
        if (hwm > capacity) fail(`MeshOps ${region} high-water mark ${hwm} exceeds capacity ${capacity}`);
      }

      // A wrapper held across clearToolingMemory() aliases reclaimed storage.
      // Using it must reject (null + STALE_WRAPPER), never abort the module —
      // the JS layer serializes tasks against exactly this ordering slip, and a
      // trap here would take the whole page down instead of one call.
      const stale = MeshOps.fromSolidName(solidName);
      if (!stale) {
        fail(`fromSolidName("${solidName}") returned null before the wipe`);
      } else {
        MeshOps.clearToolingMemory();
        for (const [what, result] of
             [['getVertices', stale.getVertices()], ['getFaces', stale.getFaces()],
              ['classifyFaces', stale.classifyFaces()], ['dual', stale.dual()]]) {
          if (result) fail(`stale wrapper ${what}() should return null`);
          if (MeshOps.getLastResult() !== MR.STALE_WRAPPER) {
            fail(`stale wrapper ${what}() did not report STALE_WRAPPER`);
          }
        }
        stale.delete();
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
      console.log(`  MeshOps: ${solidName} + dual, truncate, relax(+clamp), hankin reject, classifyFaces, clearToolingMemory, MeshOpResult reasons OK`);

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
  // These free functions and PaletteOps recipe compilation back the JS tools but are
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

  // PaletteOps V4: 256*3 sRGB bytes and engine-authored diagnostics.
  {
    const po = new Module.PaletteOps();
    try {
      const recipe = {
        schemaVersion: 4, keyCount: 3, input: { offset: 0, span: 1 },
        domain: 0, easing: 1, colorPath: 0,
        hue: {
          mode: 0, harmony: 5, direction: 0, baseTurns: 0,
          spreadTurns: 0.07, sweepTurns: 1, customTurns: [0, 0, 0, 0, 0, 0],
        },
        lightness: { curve: 1, center: 0.52, range: 0.72, custom: [0, 0, 0, 0, 0, 0] },
        chroma: {
          curve: 0, basis: 0, center: 0.62, range: 0,
          headroom: 0.94, custom: [0, 0, 0, 0, 0, 0],
        },
        hueTorsion: 0, falloffStart: 0.9,
      };
      const result = po.inspectV4(recipe);
      const lut = Uint8Array.from(result.lut ?? []);
      if (result.status.code !== 0) {
        fail(`compileAndBakeV4 returned status ${result.status.code}`);
      }
      if (!lut || lut.length !== 256 * 3) {
        fail(`compileAndBakeV4 LUT length ${lut && lut.length}, expected ${256 * 3}`);
      } else if (lut[0] === lut[765] && lut[1] === lut[766] && lut[2] === lut[767]) {
        fail(`compileAndBakeV4 gradient is flat end-to-end: [${lut[0]},${lut[1]},${lut[2]}]`);
      }
      if (!result.diagnostics || result.diagnostics.length !== 256 * 6)
        fail(`inspectV4 diagnostics length ${result.diagnostics?.length}, expected ${256 * 6}`);

      recipe.input = { offset: 0.2, span: 0.4 };
      const windowed = po.compileAndBakeV4(recipe);
      const windowedLut = Uint8Array.from(windowed.lut ?? []);
      const sourceEndpoints = [51, 153];
      for (const [destination, source] of [[0, sourceEndpoints[0]], [255, sourceEndpoints[1]]]) {
        for (let channel = 0; channel < 3; ++channel) {
          if (windowedLut[destination * 3 + channel] !== lut[source * 3 + channel])
            fail(`palette input window endpoint ${destination}, channel ${channel} did not map to source ${source}`);
        }
      }
      if (Math.abs(windowed.canonicalRecipe.input.offset - 0.2) > 1e-6 ||
          Math.abs(windowed.canonicalRecipe.input.span - 0.4) > 1e-6)
        fail(`palette input window was not preserved by canonicalization`);
    } finally {
      po.delete();
    }
  }
  console.log('  color/palette/geometry: transfer, interp, OKLab, HSV, procedural, lissajous, mobius, palette V4 OK');

  disarmDynamicExecProbe();
  if (dynamicExecCalls > 0) {
    fail(`the glue generated code ${dynamicExecCalls} time(s) — it needs ` +
      `'unsafe-eval', which the tool pages' CSP does not grant. Check the ` +
      `-sDYNAMIC_EXECUTION=0 / -sEMBIND_AOT=1 link options. First site:\n${firstDynamicExec}`);
  } else {
    console.log('\nCSP: the glue generated no code across the whole run (no eval / new Function)');
  }

  if (failures > 0) {
    console.error(`\nwasm_smoke: ${failures} failure(s)`);
    process.exitCode = 1;
    return;
  }
  console.log('\nwasm_smoke: OK');
}

await main();
