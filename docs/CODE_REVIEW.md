# Holosphere + Daydream — Code Quality Review

**Date:** 2026-06-04
**Scope:** ~50k lines C++ (engine, firmware, effects, WASM bridge, tests) across `Holosphere`, plus ~2.4k lines of vanilla-JS web simulator in `daydream`. Reviewed by 8 parallel component audits, each reading source directly and citing `file:line`. Generated tables (`reaction_graph.cpp`, `color_luts.h`) and third-party (`FastNoiseLite.h`) were judged on integration, not line-by-line.

## Overall grade: **A− / B+**

This is an unusually sophisticated and well-engineered codebase for what is effectively a solo art project. The architecture is principled (compile-time `<W,H>` specialization, 16-bit linear color, arena allocation, variadic filter pipeline), the performance engineering is genuinely expert-level (LUT-driven trig, fast-reject heuristics, zero-copy WASM readback, ISR double-buffering), and the documentation is exceptional. It is held back from an unqualified A by a **recurring class of silent out-of-bounds-on-overflow hazards** and a handful of concrete correctness bugs in the rasterizer sink, a few effects, and the WASM memory-view contract.

## Dimension grades (aggregated across all components)

| Dimension | Grade | Justification |
|---|---|---|
| **Architecture & Design** | **A−** | Compile-time resolution, clean platform abstraction, composable variadic filter pipeline, explicit-`Arena&` purity, pub/sub web state. Consistent and orthogonal across the stack. |
| **Readability & Maintainability** | **B+** | Strong Doxygen, named tolerances, intent-revealing comments. Dragged down by copy-paste (palette banks, MeshOps wrappers, magic numbers) and a few 270-line monolith functions. |
| **Correctness & Robustness** | **B / B−** | Defensive at poles/wrap/NaN in the math core; but real hazards: stripped-assert overflow → silent OOB, Voronoi OOB, div-by-zero, WASM view detachment. |
| **Memory Safety & Resource Mgmt** | **B+ / B** | Disciplined arena ownership, `Persist`/`ScratchScope` RAII, debug use-after-free tracking, Three.js GPU disposal. Undercut by unchecked `Arena::allocate` derefs and stripped capacity guards on device. |
| **Performance & Efficiency** | **A−** | Best dimension. Microsecond-budget awareness everywhere; LUTs, fast-reject culling, non-blocking DMA, instanced/zero-alloc render loops. |
| **Testing** | **B** | Excellent native harness with real algebraic invariants (Euler characteristic, partition-of-unity, error-bounded epsilons). But effect tests are near-tautological, 3 effects are quarantined rather than fixed, and integration/hardware/WASM layers are untested. |
| **Documentation** | **A−** | README is outstanding; code comments quantify error bounds and rationale. Bridge API lacks a JS-facing contract doc. |

### Per-component snapshot

| Component | Arch | Correctness | Perf | Notes |
|---|---|---|---|---|
| Core Math & Color | A− | B | A | `3dmath.h` is the strongest file in the repo |
| Rendering Pipeline | A | B | A | Sink wrap-robustness is the one real hazard |
| Animation & Mesh | A− | B | A− | Unchecked `allocate` derefs; topology assumes manifold |
| Memory/Platform/HW | A− | **B** | A− | ISR `effect_` lifetime at teardown; relaxed DMA-flag atomics |
| Effects (×27) | A− | B− | A− | Modern vs legacy quality gap; raw-loop effects carry the bugs |
| Test Suite | A− | A− | — | Strong harness, shallow on effects/integration |
| WASM Bridge & Build | A− | B− | A | View-detachment contract unenforced |
| Daydream Simulator | A− | B+ | A | Detached-buffer guard done correctly; minor worker/clip gaps |

---

## Failure philosophy (the standard these findings are measured against)

This codebase runs on hardware in a **fail-fast** posture, and the recommendations below assume it:

- **Invariant violations must hard-fault, not degrade.** An arena over-allocation, container capacity overflow, OOM, or out-of-bounds index is a logic/sizing bug. There is no valid recovery — a truncated or aliased allocation just relocates the corruption into the rendered output. A corrupted/garbage effect is a *worse* outcome than a crash, because it's silent and ships. These should crash at the violation site so they're caught on the bench.
- **The crash must survive the release build.** The hardware runs an optimized build (`NDEBUG` defined) for true performance — debug builds are not an option. `assert()` is therefore the wrong tool: it's stripped by `NDEBUG` and pulls in newlib's `__assert_func → fprintf` (the documented reason it was stripped at `platform.h:21`). The right tool is a lightweight always-on check that survives `NDEBUG` and pulls in no stdio:
  ```cpp
  #define HS_CHECK(c) do { if (!(c)) __builtin_trap(); } while (0)   // or a tiny panic/LED handler then trap
  ```
  Placed on **cold paths only** (allocation, container growth — never the per-pixel loop), it's a single predicted-not-taken branch with zero hot-loop cost.
- **Bounded/soft handling is correct only for genuine *transient* conditions** that are expected to occur and are not logic bugs — e.g. a DMA overrun dropping a frame, a missed sync pulse, last-good fallback on a sensor read. Distinguish "should never happen" (trap) from "happens occasionally by design" (degrade).
- **The host/sim/test build keeps `assert` on** — host-side violations are free to catch before they reach the device.

The current code mostly does the *opposite* of fail-fast on the invariant paths: `NDEBUG` strips the guards, so over-allocation silently writes OOB instead of crashing. The P0/P1 items below convert those silent-UB sites to `HS_CHECK` traps.

---

## Prioritized fix list

### P0 — Correctness hazards (can corrupt memory, crash a live show, or render wrong output)

1. **✅ FIXED (2026-06-04).** Added an always-on `HS_CHECK(cond)` trap (`platform.h` — `__builtin_trap()` on failure, survives `NDEBUG`, no stdio) and applied it on the cold paths: `ArenaVector::push_back/append_bulk/emplace_back` capacity guards (`memory.h`) and the `_steps_cache` guard (`plot.h`, now checked before the push) became `HS_CHECK`; `Arena::allocate` now `HS_CHECK`-traps on OOM internally **instead of returning `nullptr`**, which transitively covers the unchecked-deref sites (`conway.h`, `mesh.h`, `hankin.h`) with no per-site edits since `allocate` can no longer return null. Hot-path bounds asserts (`operator[]`, `back`) were left as debug-only `assert` to keep the per-pixel loop cost-free. `test_arena_oom_returns_null` was replaced by `test_arena_fills_to_capacity` to match the new trap-on-OOM contract. NOTE: this converts the 3 quarantined effects' *silent* fixed-capacity overflow at 288×144 into a *visible* trap, but right-sizing their capacities to actually run is tracked separately as **P2 #9**. Original finding below.

   ~~`platform.h:21` defines `NDEBUG` on device — the in-code comment shows the intent was code-size/dependency (`// Strip assert() to avoid linking newlib's __assert_func → fprintf`), **not** a deliberate fail-soft policy. The effect is that the `assert`-based capacity guards in `ArenaVector::push_back/emplace_back/append_bulk` (`memory.h:250,261,269`), `plot.h:189-191`, and the unchecked `Arena::allocate` derefs in `conway.h:206-212`, `mesh.h:135`, `hankin.h:90-95` all vanish in the performance build. On overflow/OOM the device does **neither** a crash **nor** a fallback — it silently writes past arena capacity. This is also the **root cause of the 3 quarantined effects** (SplineFlow/TestShapes/Thrusters abort at 288×144). **Fix (fail-fast, full performance, no debug build): add a lightweight always-on fault check that survives `NDEBUG` and pulls in no stdio**, e.g. `#define HS_CHECK(c) do { if(!(c)) __builtin_trap(); } while(0)` (or a tiny panic handler). Use it in place of these stripped `assert`s on the **cold paths** only (container growth, arena OOM null) — a single predicted-not-taken branch, zero hot-loop cost — so an invariant violation hard-faults immediately during testing instead of corrupting the arena. Reserve bounded/soft handling for genuine *transient* runtime conditions (DMA overrun, dropped frame), not invariant violations.~~

2. ~~**Voronoi OOB + div-by-zero** (`effects/Voronoi.h:71,123`): `sites_buffer[bestIdx]` with `bestIdx==-1` when empty; `i/(num_sites-1)` divides by zero at `num_sites==1`. Resolve by intent: if 0/1 sites is a **legitimate parameter state**, guard it (clamp/skip the degenerate term) — that's valid input, not a bug; if `num_sites` can never legitimately reach 0/1 here, `HS_CHECK(num_sites >= 2)` so a bad caller traps instead of reading OOB. Same call: GS reaction-diffusion `wb/tw` NaN (`GSReactionDiffusion.h:128`) — an empty kernel is a genuine runtime case, so guard `tw>0` (BZ already does); the NaN silently poisoning `palette.get` is the masking to avoid.~~ — **✅ FIXED (2026-06-04).** Chose the **guard** path for Voronoi (`num_sites` is not a user-registered param and the code already guards the `t` computation against `num_sites<=1`, so degeneracy is a legitimate-but-degenerate state, not a caller bug): bail before the pixel loop when `sites_buffer` is empty (no `bestIdx==-1` deref), and reuse the existing `> 1 ? n-1 : 1` span guard for the `y` term. GS: guard `if (tw <= 0.0001f) return 0.0f;` — mirrors BZ's `sample_kernel`, and the `0.0f` falls through the shader's `b < 0.05f` cull to transparent black.

3. ~~**WASM detached-buffer hazard** (`wasm.cpp:101,169-172,207-219`). With `ALLOW_MEMORY_GROWTH=1`, any allocation can detach the ArrayBuffer behind `getPixels()`/`getParamValues()` views. Pre-size `pixelBuffer` to `MAX_W*MAX_H*3` once and never `resize()`; re-fetch the view each frame; document the contract. (Daydream already guards this correctly at `daydream.js:357` — mirror that expectation in the bridge.)~~ — **✅ FIXED (2026-06-04).** `pixelBuffer` is pre-sized to `MAX_W*MAX_H*3` and `paramValues` pre-reserved to `MAX_PARAMS` (32) once in the `HolosphereEngine` ctor; `setResolution` no longer `resize()`s; `getPixels()`/`getBufferLength()` now return a view/length over the active-resolution prefix of the stable buffer. The steady-state render/sync path performs no allocation, so the only detachment source is unrelated heap growth — already guarded by `daydream.js::refreshPixelView`. The memory-view contract is documented on `getPixels()`/`getParamValues()`.

### P1 — Robustness, concurrency & lifetime

4. ~~**ISR reads `effect_` with no publish/lifetime barrier** (`pov_single.h:78-84`, `pov_segmented.h`). A late-pending IntervalTimer ISR can dereference a freed effect at teardown. Null-publish `effect_` under disabled interrupts before `delete`.~~ — **❌ DISMISSED (2026-06-04), not a bug as described.** The lifetime barrier already exists: `run()` calls `timer.end()` (`pov_single.h:108`) / `detachInterrupt()` (`pov_segmented.h:261-262`) to stop the ISR *before* control returns to `show()` and runs `delete effect_`. Any genuinely pending interrupt is serviced the instant interrupts resume — right after `end()`/`detach()` and well before the `delete` — so it reads a still-valid `effect_`, not freed memory. The proposed fix is also incoherent for this code: **neither ISR null-checks `effect_`** (`pov_single.h:124-153`, `pov_segmented.h:282-317` deref it unconditionally), so null-publishing before `delete` would convert a hypothetical use-after-free into a *guaranteed* null-deref on any stray ISR. The only real kernel — the delete-before-null ordering and an unguarded ISR-shared raw pointer — is defense-in-depth, not a live hazard, and any real fix would need an ISR null-check **plus** reordering (not the doc's one-liner); deferred as not worth the per-column ISR branch cost.

5. ~~**Relaxed atomics on the DMA completion flag** (`dma_led.h:270-295`). Use `release` in the ISR / `acquire` on the consumer instead of `relaxed` to guarantee buffer-write→flag visibility.~~ — **❌ DISMISSED (2026-06-04), not a bug.** Teensy 4.x is a **single-core** Cortex-M7: the DMA-completion ISR and the consumer thread are the *same* observer (the ISR preempts), so acquire/release vs. relaxed cannot change what they see — relaxed is correct (the same single-core rationale the code already documents for `instance_` at `dma_led.h:298-306`). The finding's stated goal, "buffer-write→flag visibility", points at the wrong mechanism: the buffer is consumed by the **DMA engine**, a separate bus master that is not a C++ observer, so the atomic's memory order has no effect on it — that coherence is provided by `arm_dcache_flush_delete()` (cache flush + DSB) before `dma_.enable()` (`dma_led.h:127/173`). Surrounding compiler reordering is already blocked by those cache-op/register-write barriers. acquire/release would only add `DMB` cost for no correctness gain on this target. Added a clarifying comment at `transferComplete_` to pre-empt re-flagging.

6. ~~**WASM `setResolution` leaves the engine dead** on an unsupported size (`wasm.cpp:87-127`): nulls `currentEffect` and renders blank with no error to JS. Centralize the supported `(W,H)` dispatch table (currently duplicated 3×) and surface failure.~~ — **✅ FIXED (2026-06-04).** Added a single `HS_WASM_RESOLUTIONS` X-macro table (`wasm.cpp`); `setResolution`/`setEffect`/`getEffectSizes` now dispatch through it so the supported set can't drift. `setResolution` rejects unsupported sizes (instead of switching to them and nulling `currentEffect`), keeps the prior valid resolution/effect alive, and returns `bool`; `daydream.js::applyResolution` checks that bool and bails (keeping the current resolution) rather than driving the engine blank.

7. ~~**Geometry edge cases:** `AABB::intersectRay` div-by-zero on axis-aligned rays (`spatial.h:64-97`); `KDTree::nearest` unclamped `k>MAX_K` silently drops results (`spatial.h:152-209`); `OrientationTrail::get` has contradictory oldest/newest ordering docs (`animation.h:230-233`) that can reverse motion-blur direction.~~ — **✅ FIXED (2026-06-04).** `intersectRay` rewritten as a per-axis slab test that guards a zero direction component (parallel ray ⇒ hit only if the origin is inside the slab), eliminating the `0/0 = NaN` on a grazing on-face ray; added `test_aabb_ray_parallel_grazing` to lock it in. `KDTree::nearest` now `HS_CHECK(k <= MAX_K)` so an over-large `k` traps instead of silently truncating (asking for more neighbors than points exist still legitimately caps below). `OrientationTrail::get` docs corrected to match the implementation — **0 is the oldest** snapshot, `length()-1` the newest (the inline comment was right; the `@param` was wrong).

8. **Worker clip not re-applied on resolution change** (`segment_worker.js:79-88`) — segments render with a stale clip rectangle until the next effect switch. Call `applyClip()` after `engine.setResolution`.

### P2 — Quality, maintainability & test depth

9. **Strengthen effect tests** (`test_effects.h:84-87`): `HS_EXPECT_TRUE(acc+1>0)` is tautological. Assert non-trivial output (not-all-zero / finite bounds) or golden-hash buffers; inject a mockable clock for determinism; root-cause the 3 quarantined effects rather than excluding them.

10. **De-duplicate:** copy-pasted `BakedPaletteBank` block (`HankinSolids.h` ↔ `IslamicStars.h`), the 15× MeshOps wrappers (`wasm.cpp:280-596`), repeated magic numbers (`62`ms in daydream, `D_AVG=0.04044f`, Hopf tumble constants), and the dual URL-sync layers (`state.js` vs `gui.js`).

11. **Unclamped lerps** (`color.h:366-368` `lerp8`, `color.h:584` `Gradient::get`) can wrap to garbage for `t∉[0,1]`; the `Face` SDF hardcodes `W=288` in its dist-LUT threshold (`sdf.h:1136`), breaking multi-resolution; decompose the 270-line `Face` ctor.

12. **Cleanup:** per-frame `Serial.print` in the POV hot loop (`pov_single.h:105-106`); dead/commented code in the `.ino` files and `effects_legacy.h`; `THREE_VERSION` mismatch (`tools/shared.js:105` vs `vendor-importmap.js:16`).

13. **Build fragility:** the install step writes into the sibling `../daydream` checkout (`CMakeLists.txt:62-66`) with no `EXISTS`/`OPTIONAL` guards — guard it or expose a `DAYDREAM_DIR` cache var.

**Notable coverage gaps to close over time:** `canvas.h`, `scan.h` rasterizer, `transformers.h`, `generators.h`, the hardware drivers, and the WASM bridge are entirely untested.

---

## Component detail

### Core Math & Color (`3dmath.h`, `color.h`, `palettes.h`, `easing.h`, `waves.h`, `util.h`, `geometry.h`, `concepts.h`, `rotate.h`, `color_luts.h`)

Strongest cluster. `3dmath.h` features branch/div-free quaternion rotation with documented op counts, fast trig/acos/atan2 with **measured** error bounds, and singularity sentinels throughout. Color path correctly blends in 16-bit linear light with an exact `/65535` reduction and an ARM SMLAD fast path; OKLCH interpolation uses shortest-arc hue.

- **High** | `color.h:366-368` — `lerp8` casts to `uint8_t` with no clamp; `t∉[0,1]` (producible by Breathe/Ripple modifiers) overflows to garbage colors.
- **Medium** | `color.h:584`, `geometry.h:281` — `Gradient::get` / `vector_to_pixel` float-`wrap` TOLERANCE-snap can drop the last column or snap valid small x to 0.
- **Medium** | `3dmath.h:859-872` — `fast_sinf` range reduction can leave x marginally out of `[0,π]`, sign-flipping a sample feeding slerp/make_rotation. Clamp post-reduction.
- **Low** | `color.h:1182-1192` — `BakedPalette::get` derefs `lut_` with no null guard before `bake()`.

### Rendering Pipeline (`filter.h`, `sdf.h`, `scan.h`, `plot.h`, `canvas.h`, `styles.h`, `transformers.h`)

Clean compile-time `Pipeline<W,H,Filters...>` recursion with automatic 2D↔3D domain bridging; SDF CSG composes orthogonally; pervasive fast-reject culling and pole guards.

- **Medium** | `filter.h:95-104` — the float `plot` sink uses `fast_wrap` (corrects only one period); a future filter emitting `x≥2W` writes OOB via `cv(xi,y)` with the assert stripped. Use modulo `wrap()` in the sink.
- **Low** | `sdf.h:1136` — `Face` dist-LUT threshold hardcodes `W=288`, mis-scaling on other resolutions.
- Nits: duplicate `#include "memory.h"` (`filter.h:23`); `DistortedRing` uses `cosf` where `Ring` uses `fast_cosf`.

### Animation & Mesh/Geometry (`animation.h`, `mesh.h`, `conway.h`, `hankin.h`, `solids.h`, `spatial.h`, `generators.h`, `presets.h`, `static_circular_buffer.h`, `reaction_graph.*`)

Type-erased inline animation storage, correct half-edge construction, KNN graph as a `PROGMEM` table with on-the-fly node reconstruction.

- ~~**High** | `Arena::allocate` returns `nullptr` on OOM but Conway/Hankin/HalfEdge scratch allocations deref unconditionally (`conway.h:206-212`, `mesh.h:135`, `hankin.h:90-95`). Per fail-fast, `allocate` should `HS_CHECK` non-null internally (over-allocation = bug → trap) rather than returning null for callers to ignore.~~ **✅ FIXED (2026-06-04)** — `allocate` HS_CHECK-traps on OOM (see P0 #1).
- ~~**High** | `ArenaVector` capacity guards are stripped `assert`s on device (`memory.h:250,261,269`) → silent OOB. Replace with `HS_CHECK` so a capacity overrun traps in the release build instead of corrupting the arena.~~ **✅ FIXED (2026-06-04)** — growth guards are now `HS_CHECK` (see P0 #1).
- **Medium (partially fixed)** | ~~`spatial.h:64-97` ray-AABB div-by-zero~~ **✅ FIXED**; ~~`spatial.h:152-209` unclamped `k` (a `k>MAX_K` caller is a bug → `HS_CHECK`, not silent drop)~~ **✅ FIXED**; ~~`animation.h:230-233` contradictory trail-ordering docs~~ **✅ FIXED** (all 2026-06-04, see P1 #7); `conway.h:78-104` `vertex_orbit` assumes manifold closure (add a hard `count<=I` loop bound independent of asserts) — **still open**.
- **Low** | `animation.h:1622` uses `Serial.println` directly (breaks platform abstraction, won't compile on WASM); `static_circular_buffer.h:124-168` already does `abort()`-on-misuse — that's the *right* fail-fast behavior, but its header comment ("No asserts — prevents hard crashes in live shows") is now misleading and should be corrected; route the `abort()` through the same `HS_CHECK`/trap path for consistency and to avoid stdio.

### Memory, Platform & Hardware (`memory.*`, `platform.h`, `led.h`, `dma_led.h`, `pov_single.h`, `pov_segmented.h`, `*.ino`)

Non-blocking DMA with a composite black-frame trick; minimal interrupt-disabled critical sections; branchless per-segment ISR resolved at boot.

- **Medium** | `memory.cpp:28-42` — `configure_arenas` silently proportionally **scales down** an over-subscribed partition request. Under fail-fast this masks a config/sizing bug: the effect asked for a layout that doesn't fit and gets a quietly smaller one, which then over-runs later. Prefer `HS_CHECK(persistent + a + b <= GLOBAL_ARENA_SIZE)` so a bad partition request traps at `init()` instead of degrading. (Previously praised as "graceful" — that was the old fail-soft framing.)
- ~~**High** | `pov_single.h:78-84,124-153` — ISR reads `effect_` with no publish/lifetime barrier.~~ **❌ DISMISSED (2026-06-04)** — `timer.end()`/`detachInterrupt()` already stop the ISR before `delete`; proposed null-publish fix is incoherent (ISR has no null-check). See P1 #4.
- ~~**Medium** | `dma_led.h:270-295` — `relaxed` atomics on the cross-context completion flag.~~ **❌ DISMISSED (2026-06-04)** — single-core Cortex-M7: ISR and consumer are the same observer, so relaxed is correct; buffer↔DMA coherence is the cache flush, not the atomic. See P1 #5.
- **Medium** | `pov_single.h:105-106` — per-frame `Serial.print` in the hot loop; dead `delay(125)`.

### Effects Library (27 effects + registry/glue)

Self-registering `REGISTER_EFFECT` factory; per-effect documented memory budgets; disciplined `Persist`/`ScratchScope` use. Clear quality gap between modern engine-based effects (BZ, Hopf, Raymarch — exemplary) and raw-pixel legacy effects (which carry the NaN/OOB bugs). `BZReactionDiffusion.h` is the reference template.

- ~~**High** | `Voronoi.h:71,123` OOB + div-by-zero.~~ **✅ FIXED (2026-06-04)** — empty-buffer bail + `num_sites==1` span guard.
- ~~**Medium** | `GSReactionDiffusion.h:128` NaN from `wb/tw` (BZ guards this, GS does not — same algorithm, divergent robustness).~~ **✅ FIXED (2026-06-04)** — `tw <= 0.0001f` guard mirroring BZ.
- **Medium** | SplineFlow/TestShapes/Thrusters quarantined for fixed-capacity overrun at 288×144.
- Recurring: copy-pasted palette-bank block, duplicated magic constants, two coexisting rendering idioms (engine rasterizers vs raw double-for loops — the latter carry the bugs).

### Test Suite (20 test headers + harness)

Clean header-only harness with delta-based per-module tallies and typed domain macros (`HS_EXPECT_VEC/QUAT/COMPLEX`). Genuinely strong invariants: Euler V−E+F=2, winding/edge reciprocity, partition-of-unity for filters, error-bounded float epsilons.

- **Medium** | `test_effects.h:84-87` — tautological postcondition (`acc+1>0`).
- **Medium** | `test_effects.h:119-134` — suite is green by excluding 3 known-failing effects.
- **Medium** | no determinism seam for `hs::millis` — effects can't be regression-tested.
- Untested: `canvas.h`, `scan.h` rasterizer, `transformers.h`, `generators.h`, platform/hardware, WASM bridge.

### WASM Bridge & Build (`wasm.cpp`, `CMakeLists.txt`, presets, toolchain, build scripts)

Engine kept pristine (bindings only wrap); templated compile-time factory dispatch; dedicated tooling arena; correct `-O3 -ffast-math -msimd128 -flto` flags; ccache + toolchain auto-detection.

- ~~**High** | `wasm.cpp:101,169-172,207-219` — `typed_memory_view` detachment under memory growth.~~ **✅ FIXED (2026-06-04)** — stable pre-sized readback buffers, no per-frame allocation, contract documented.
- ~~**High** | `wasm.cpp:87-127` — `setResolution` on an unsupported size leaves the engine null/dead with no error to JS.~~ **✅ FIXED (2026-06-04)** — rejects unsupported sizes, returns bool, single X-macro dispatch table. See P1 #6.
- **Medium** | supported resolutions hardcoded in 3 places; install into sibling `../daydream` is unguarded.
- **Low** | 15× MeshOps wrapper boilerplate; 8 KB stack with no release `STACK_OVERFLOW_CHECK`.

### Daydream Web Simulator (vanilla ES-module JS + Three.js + Web Workers + WASM)

Clean pub/sub state, worker isolation, near-zero per-frame allocation (persistent adapter, pooled labels, reused temp vectors), correct GPU disposal on resolution change, and — notably — a **correct** detached-buffer guard (`daydream.js:357`). Robust browser-compat (codec fallback chain, importmap probe, File System Access fallback).

- **Medium** | `segment_worker.js:79-88` — `setResolution` doesn't re-apply clip → stale clip rectangle.
- **Medium** | `daydream.js:631,672` — `applyResolution(true)` called twice at startup (once with null engine).
- **Low** | new/resized segment workers don't receive current tuned param values; dual URL-sync layers can race; `62`ms magic number in 3+ places; `THREE_VERSION` mismatch.

---

## Bottom line

The engine reflects strong systems-design instincts and real domain expertise — the architecture, performance work, and documentation would not be out of place in a professional embedded-graphics codebase. The dominant risk is a single systemic pattern: capacity/OOM failures that are correctly *detected* in debug but silently become out-of-bounds writes in the performance build, because `NDEBUG` strips the guarding `assert`s. The intended posture is **fail-fast**: an arena over-allocation, capacity overflow, or OOM is a logic bug that should hard-fault immediately during testing — a corrupted/aliased arena producing garbage effects is the worst outcome. Addressing **P0 #1** with a lightweight always-on `HS_CHECK`→trap (compiled into the release build, no stdio dependency, cold-path only) closes the largest correctness gap and eliminates the quarantined-effects problem at full performance — without resorting to a debug build.
