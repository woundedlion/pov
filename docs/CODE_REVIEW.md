# Holosphere + daydream — Code Quality Review

**Date:** 2026-07-10
**Repositories:** `Holosphere` (C++ engine + Teensy firmware + WASM target) and `daydream` (web simulator)
**Method:** Orchestrated multi-agent review. Eighteen independent expert reviewers (one per component group) read every in-scope file in full, assessed all software-quality dimensions against the architecture described in the README, and validated each finding against the cited source before inclusion (real defect + concrete failure scenario + correct minimal fix). Speculative, cosmetic-only, already-handled, and by-design items were dropped at the agent level.
**Out of scope (per instruction):** `core/engine/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`, `core/math/rotate.h`, and the KiCad/PCB generation under `hardware/phantasm/`.

---

## Executive Summary

**Overall grade: A− (high, bordering A).**

This is exceptionally mature, heavily-audited code. Across ~140 source files in two languages, the review surfaced **zero Critical and zero High severity defects**, **2 Medium**, and **43 Low** — and the great majority of the Low items are consistency, observability, or documentation polish rather than latent bugs. Every component graded **A or A−**. The engineering discipline is uniformly high: airtight arena/memory-lifetime reasoning, compile-time-enforced invariants (filter ordering, arena budgets, buffer bounds, roster alignment), a fail-fast philosophy applied consistently on cold paths only, host/device determinism treated as a first-class contract, and a genuinely deep test apparatus (native death/determinism/parity harnesses plus 328 passing Node tests with real WASM-parity assertions).

The two Medium findings are a single forgotten rendering flag that produces a real visual seam under segmented rendering (`PetalFlow`), and a suppressed keyboard focus ring in the standalone designer tools (a WCAG 2.4.7 accessibility defect). Both are one- or few-line fixes.

The dominant theme is that the codebase is at a **polish ceiling**: what remains after prior review passes is naming drift, `assert`-vs-`HS_CHECK` consistency, missing-but-non-triggered guards, telemetry granularity, and accessibility gaps in the ancillary tool pages. There is no structural debt, no architectural rework indicated, and no correctness hot-spot.

---

## Grade Summary by Dimension (project-wide)

| Dimension | Grade | Rationale |
|---|---|---|
| Correctness | **A** | No Critical/High defects survived validation. One Medium visual-seam bug and a handful of latent-only Lows against otherwise rigorous numeric/topology/concurrency reasoning. |
| Memory safety | **A** | Arena bump/rewind invariants are wrap-proof and exact; fixed-capacity containers, generation-tracked use-after-free detection, no heap on the render path. A couple of `assert`-vs-`HS_CHECK` and enforcement-gap nits. |
| Architecture & design elegance | **A** | Standout: compile-time filter ordering via `static_assert`, the partitioned-arena model with explicit `(target, temp)` plumbing, the pure-core/device-shell split in hardware, and a codebase deliberately architected for host-side testability. |
| Interface expressiveness | **A−** | Zero-alloc type-erased callables, deleted-rvalue borrow enforcement, concept-gated CSG, and fluent builders are excellent; docked for `Effect` API naming drift, a projection-specific `operator/`, and a few footguns. |
| Readability | **A−** | Precise naming and exceptional doxygen/JSDoc; a few very long numeric routines and two monolithic daydream DOM modules demand real effort. |
| Error handling | **A−** | Disciplined fail-fast on cold seams, soft-degrade on designed overflow; docked for `assert`-vs-`HS_CHECK` inconsistency and unguarded embind calls on the JS runtime switch path. |
| Performance | **A** | Measured structural optima, baked LUTs, zero-copy readback, branchless ISR packing; only setup-path O(n²) and minor per-frame JS allocation churn remain. |
| Testability | **A** | Native death-harness, cross-run determinism, WASM-parity + cross-implementation firmware checks, Node logic extraction; monolithic `daydream.js`/`driver.js` and a WASM-smoke coverage gap keep it off A+. |
| Consistency | **A−** | Strong shared idioms; naming drift, epsilon divergence, and the one PetalFlow flag omission are the wrinkles. |
| Documentation | **A** | Among the best-documented codebases of its kind; load-bearing invariants captured at the point of use. Two inaccurate comments. |
| Accessibility (daydream) | **B** | `index.html` is strong (canvas `role="img"`, tested keyboard nav); the standalone tool pages suppress the focus ring and omit ARIA text alternatives. |

---

## Grade Summary by Component

| Component | Grade | Findings |
|---|---|---|
| core/engine — platform/memory/registry | A | 2 Low |
| core/engine — generators/transformers/styles/presets/reaction_graph | A | 1 Low |
| core/math | A− | 2 Low |
| core/mesh — mesh/mesh_classes/spatial | A− | 5 Low |
| core/mesh — conway/hankin/solids | A | 2 Low |
| core/color | A | 2 Low |
| core/render — canvas/filter/shading/led | A− | 3 Low |
| core/render — sdf/scan/plot | A− | 3 Low |
| core/animation | A | 4 Low |
| hardware — drivers/sync | A | 2 Low |
| effects — group A (RD/particle/field) | A | 2 Low |
| effects — group B (mesh/segue/mobius) | A | 1 Low |
| effects — group C (petal/ray/voronoi/…) | A− | 1 Medium |
| targets / scripts / tests | A | 2 Low |
| daydream — core JS (app/driver/state) | A− | 3 Low |
| daydream — core JS (workers/recorder/sequencing) | A | 2 Low |
| daydream — UI + geometry tools | A− | 3 Low |
| daydream — tests + HTML | A | 1 Medium, 2 Low |

---

## Notable Strengths

- **Compile-time invariant enforcement everywhere it matters.** Filter-pipeline ordering (terminal-last, world-before-screen, unit-input) is a `static_assert` chain; CSG span buffers are sized by a recursive `sdf_max_spans<>` that fails to compile on overflow; per-effect persistent arena footprints are `static_assert`-budgeted against the device partition; the effect roster is triple-pinned across registry, native smoke, and CMake. Whole classes of misuse are build errors rather than runtime faults.
- **Memory model that is both fast and provably safe.** The 330 KiB partitioned arena with explicit `(Arena& target, Arena& temp)` plumbing makes lifetime visible at every call site; bump/rewind bounds are written in wrap-proof subtractive form; `ArenaVector`/`ArenaSpan` carry generation + rebind tracking that catches use-after-free and reallocation-under-borrow in debug.
- **Host/device determinism as a contract.** Shared `Pcg32(1337)`, FastLED integer-mock parity, narrowed `millis()`/`micros()` wrap semantics, and a `__FINITE_MATH_ONLY__` `#error` guarding the NaN→hi clamp keep the WASM simulator and Teensy firmware bit-comparable.
- **A test apparatus that verifies its own safety nets.** The native suite spawns subprocesses to confirm each `HS_CHECK` actually traps (`SIGILL`), re-renders every effect under a fixed clock to diff for determinism, and pins the WASM param-marshaling streams index-aligned. daydream's 328 Node tests include real WASM-vs-JS parity assertions each backed by an independent absolute golden, plus a cross-implementation check of the JS tiling against a hand port of the C++ firmware map.
- **The Phantasm 1-wire sync core.** A position-from-time flywheel that makes masked-IRQ windows structurally unable to drop columns, an odd-only distance-2 symbol alphabet that degrades any glitch to "missed, never misclassified," and a race-closed effect-lifecycle handoff — all extracted into a pure, host-tested core.
- **Documentation that captures the landmines.** Composition polarity, tap-anchor invariants, seam-split ordering, arena contracts, and accepted numerical limits are documented at the point of use, not in a separate design doc that rots.

---

## Prioritized Remediation List

Every validated defect is listed below under its priority band, numbered in a single sequence. Each item is tagged with its component and cited `file:line`, followed by the recommended minimal fix. None require a performance trade-off.

### Priority 0 — Critical

_None._

### Priority 1 — High

_None._

### Priority 2 — Medium

1. ✅ **[effects] `effects/PetalFlow.h:47`** — Seam claim is a false positive: `AntiAlias` is a forward quintic splat with `crosses_segments = false` (inherits `Is2D`, no history), so it never gathers across a band boundary and produces no seam under segmented rendering. For PetalFlow's `Orient + AntiAlias` pipeline `any_crosses_segments` folds to `false`, so the flag is runtime-identical to the default. Set anyway for trait-derivation parity with siblings (RingShower/ShapeShifter/SplineFlow/Thrusters): a future history-bearing filter now auto-sets `full_frame` instead of silently leaving it false. *Fix applied:* base init is `Effect(W, H, {.strobe = true, .full_frame = decltype(filters)::any_crosses_segments})`.

2. ✅ **[daydream] `daydream/tools/tools.css:119-121`** — `.toggle-switch:focus { outline: none; }` strips the focus ring from the keyboard-focusable `<button role="switch">` toggles (solids.html) with no replacement, so a keyboard user tabbing the toggles gets no visible focus indicator (WCAG 2.4.7). *Fix:* replace with `.toggle-switch:focus-visible { outline: 2px solid var(--blue-500); outline-offset: 2px; }` and drop the blanket `:focus { outline: none }`.

### Priority 3 — Low

_Latent correctness / robustness (real behavior under a reachable condition):_

3. ✅ **[animation] `core/animation/timeline.h:305-317`** — `new_vals_count = global_timeline_num_events - active_cnt` assumes callbacks only *add* events; a `.then()`/timer callback that calls `Timeline::clear()` mid-step drives the count negative and corrupts subsequent steps. *Fix:* add `HS_CHECK(global_timeline_num_events >= active_cnt)` after the loop, or document that callbacks must not clear/shrink the timeline.

4. ✅ **[animation] `core/animation/params.h:685-693`** — `Ripple` ctor snapshots `peak_amplitude(params.amplitude)` then zeroes `params.amplitude`; a caller who sets `params.amplitude` *after* constructing gets a silent no-op ripple. *Fix:* take the peak as an explicit ctor argument, or document that `params.amplitude` is read at construction and must be set beforehand.

5. ✅ **[daydream] `daydream/daydream.js:461-477`** — The runtime effect/resolution-switch path calls `setEffect`/`getParameterDefinitions`/`setResolution` unguarded, unlike the initial-load path (try/catch at line 567). A WASM abort/`BindingError` (rather than the expected `false` return) propagates through `AppState.notify`, leaving state half-updated. *Fix:* wrap the `applyEffect()`/`applyResolution()` bodies in the same try/catch used on load.

6. ✅ **[daydream] `daydream/daydream.js:650-668`** — The segmented-toggle `onChange` handlers re-check `segments.active` after `await warmModules()`, but an on→off→on burst leaves two in-flight handlers that both pass the re-check and both call `segments.create()`, double-spawning the worker pool. *Fix:* capture a monotonic epoch/token before the await and compare after, or serialize create/destroy.

7. ✅ **[daydream] `daydream/recorder.js:383-386`** — Cancelling the Save-As dialog sets `aborted=true` but does not stop the already-started `MediaRecorder`; the session keeps capturing while every chunk is discarded, and stop yields no file with only a console warning. *Fix:* call `this.stop()` in the `opened.catch` `AbortError` branch so a cancelled save ends the session.

8. ✅ **[color] `core/color/palettes.h:161,167`** — `MeshPaletteBank::operator[](int i)` indexes the 5-element array with no bounds trap, against the codebase's fail-fast convention (all sibling palette accessors check). Safe only because current indices come from `shuffle_indices`. *Fix:* add `HS_CHECK(i >= 0 && i < N, ...)` to both overloads (off the per-pixel path).

9. ✅ **[engine] `core/engine/inplace_function.h:137-149`** — `ArenaVector`'s destructor-skipping contract relies on captured callables being trivially destructible, but `inplace_function`'s converting ctor only `static_assert`s copy/nothrow-move/nothrow-copy — never trivial destructibility. A nothrow-copyable callable owning an external resource could be stored and silently leaked. *Fix:* add `static_assert(std::is_trivially_destructible_v<D>, ...)` to the converting ctor (all current captures already satisfy it).

10. ✅ **[mesh] `core/mesh/mesh.h:663-665`** — `classify_faces_impl` rebinds `mesh.topology` to `persistent` only when `capacity() < F`; a second classify passing a *different* `persistent` arena while capacity already suffices silently keeps topology in the original arena. *Fix:* assert the bound arena matches `persistent`, or bind unconditionally when the arena identity differs.

11. ✅ **[mesh] `core/mesh/spatial.h:82-89`** — `nodes.bind(arena, count)` runs before `HS_CHECK(count <= MAX_POINTS)`, so an oversized point set traps as a generic arena-overflow inside `bind` rather than at the specific int16 node-link diagnostic. *Fix:* move the `MAX_POINTS` check above the `bind`.

12. ✅ **[render] `core/render/scan.h:65-83,828`** — `Union`/`SmoothUnion` set `is_solid = A::is_solid || B::is_solid`, so `process_pixel` picks the solid AA branch even when the winning (nearer) child is a stroke; e.g. `Union<Circle, Ring>` renders the ring child as a hard band instead of its stroke falloff. *Fix:* route the winner's `size`-based falloff when the winning child is non-solid, or document that mixing solidity classes under one CSG node is unsupported.

13. ❌ **[render] `core/render/sdf.h:2222-2258`** — Real but latent, rejected as disproportionate. Azimuth is monotonic along any great-circle arc that does not cross a pole, so per-vertex thetas already bound each edge's column coverage; the gap arises only for an edge grazing *just outside* a pole (endpoints near-antipodal in longitude), which no shipped mesh produces — its polar faces put a vertex *at* the pole or *enclose* it (handled by `apply_pole_containment`). The proposed `full_width` fix over-triggers on the common vertex-at-pole faces (their pole-adjacent edge extrema sit at phi≈0), flipping large polar faces to full-width column scans — a real perf cost to guard a configuration mesh generation cannot generate.

14. ✅ **[effects] `effects/Comets.h:156-161`** — `closing_domain` floors `m2·domain` but not the divisor `config.m2` itself, on a table the class doc explicitly frames as an authored extension point; a future entry with `m2 == 0` yields an infinite domain that freezes the head. *Fix:* floor `m2` or `HS_CHECK(config.m2 > 0)`.

15. ✅ **[effects] `effects/MindSplatter.h:277-284`** — `assert(p_idx < active_count)` runs *before* the defensive clamp the adjacent comment justifies as float-overshoot protection; on builds where the assert is live the clamp is dead, and on device (assert stripped) the safety net is never checked. *Fix:* clamp first, then assert on the clamped value (or drop the assert and keep the `active_count`-guarded clamp).

16. ✅ **[daydream] `daydream/tools/banner.js:35`** — `showFatalError` appends the banner only when `document.body` is truthy; a fatal error before `<body>` exists builds the element but never inserts it, and later calls' `getElementById` return null, so the banner is silently lost. *Fix:* `(document.body || document.documentElement).appendChild(el)`.

17. ❌ **[engine] `core/engine/memory.h:538-543,576-591`** — Rejected: `ArenaVector::operator[]` is on the per-pixel path (`core/render/plot.h:653` indexes `steps_cache[j]` once per plotted point, ~48k plots/frame), so promoting `assert`→`HS_CHECK` would add an always-on branch per plotted pixel — a perf regression. The `assert` is the intentional hot-access choice; the sibling `StaticCircularBuffer` is not on that path.

18. ✅ **[effects] `effects/DreamBalls.h:57-71`** — Unlike every sibling (Comets/Dynamo/FlowField/BZ/GS), `DreamBalls` has no compile-time persistent-arena budget `static_assert`; adding a preset or a larger solid would overrun the default partition as a runtime trap instead of a build error. *Fix:* `static_assert(2 * BakedPalette::LUT_SIZE * sizeof(Color4) <= PERSISTENT_BUDGET, ...)` with a note that the mesh geometry is runtime-bounded.

_Observability, diagnostics, and test coverage:_

19. ✅ **[tests] `scripts/wasm_smoke.mjs`** — The CI embind smoke harness exercises `HolosphereEngine`/`MeshOps`/the spline exports but never touches `PaletteOps.bakeLut` or the 8 color/lissajous free-function exports registered in `wasm.cpp`, so a wrong-target/transposed-arg binding there ships green through Holosphere CI (only daydream's separate parity tests cover it). *Fix:* add one call + finite-result assertion per uncovered export, mirroring the spline block.

20. ✅ **[hardware] `hardware/pov_sync.h:1272,1555`** — `telemetry_.emit_aborted` is incremented for two distinct events (mid-burst lateness self-censor vs. dropping a stale overrun beacon), but its doc describes only the first, so a rising `abrt` in field telemetry conflates two root causes. *Fix:* add a distinct `beacons_overrun_dropped` counter, or amend the doc to state it counts both.

21. ✅ **[daydream] `daydream/tools/*.html` (lissajous:39, mobius:254, solids:471, splines:75)** — The interactive 3D viewport `<canvas>` elements have no text alternative, unlike `index.html` (`role="img"` + `aria-label`); a screen reader announces a bare "canvas." *Fix:* add an `aria-label` naming each viewport.

22. ✅ **[daydream] `daydream/tools/splines.html:94,99` and `palettes.html:229,327`** — Range inputs are labeled only by an adjacent `<span>` with no programmatic association. *Fix:* give the label span an `id` and add `aria-labelledby` (or `aria-label`) to the input.

_Consistency and maintainability:_

23. ✅ **[render] `core/render/canvas.h` (397/417/427/471 vs 458/484/105)** — The public `Effect` API mixes camelCase (`updateParameter`, `getParameters`, `setAnimationsPaused`, `markReadonly`) with snake_case (`register_param`, `mark_animated`, `needs_full_frame`, `set_clip`) in one widely-derived base class. *Resolved:* documented the split at the `Effect` class doc — the camelCase names are the JS/embind boundary (bound in `wasm.cpp`), snake_case is the internal C++ API; renaming would break the WASM bridge.

24. ✅ **[render] `core/render/led.h:39-43`** — The `USE_DMA_LEDS` stubs `NoColorCorrection`/`NoTempCorrection` are empty/trivial, whereas the real guards delete copy and have non-trivial dtors, so copy-init call sites compile only in the non-DMA build and unused stubs risk `-Wunused-variable`. *Fix:* give the stubs `= default` dtors + deleted copy for identical semantics, or mark call sites `[[maybe_unused]]`.

25. ✅ **[render] `core/render/filter.h:1020 vs 1203,1192`** — Near-zero splat weight cutoffs diverge: `AntiAlias` uses `1e-8f`, `Blur` uses `1e-5f`. No artifact, but the magic numbers obscure intent. *Resolved:* commented why they differ at each site (raw bilinear coverage products vs. normalized 3x3 kernel taps) — not unified, the divergence is intentional.

26. **[mesh] `core/mesh/mesh.h:360-382`** — `compile()`'s signature advertises only `geom_arena` yet silently consumes global `scratch_arena_a`, breaking the explicit-`(target, temp)` convention that `classify_faces_by_topology` follows and hurting reentrancy/testability. *Fix:* accept an explicit scratch `Arena&`.

27. **[mesh] `core/mesh/mesh.h:503`** — `classify_faces_impl` is `always_inline` but is a ~180-line cold routine force-inlined into two `HS_COLD` wrappers, duplicating cold code with no perf benefit on the memory-tight Teensy image. *Fix:* drop `always_inline` and emit once.

28. **[mesh] `core/mesh/solids.h:1342-1348`** — `build_vertex_directions` computes per-vertex nearest-neighbor angle by brute-force O(n²) dot products (up to ~76M iterations for a large solid) despite the codebase already shipping a KD-tree for exactly this query. Setup-only/`HS_COLD`, hence Low. *Fix:* build the KD-tree once and query nearest neighbors, or comment that the O(n²) is deliberate for the bounded setup path.

29. **[mesh] `core/mesh/conway.h:660` vs `core/mesh/solids.h:279`** — The `expand` default factor is spelled two ways: `2.0f - sqrtf(2.0f)` (runtime call in a default argument) vs `2.0f - SQRT2` (constexpr). Values match but invite drift. *Fix:* hoist a shared `constexpr float EXPAND_DEFAULT_T` referenced by both.

30. **[math] `core/math/3dmath.h:647-658`** — `Complex::operator/` is not general complex division: it silently caps magnitude at `STEREO_INF` and returns `(0,0)` for 0/0. Correct/intentional for the stereographic call sites, but the operator carries no name signaling it, so any future generic use gets silently clamped quotients. *Fix:* rename to a free `project_div()` used by `mobius`/`stereo`, or add a "projection-domain only" caveat at the operator itself.

31. **[animation] `core/animation/sprites.h:169-189`** — `ParticleSystem` exposes `pool`, `active_count`, `friction`, `gravity`, `max_life`, `attractors`, `emitters` as public mutable members, unlike every other animation (all-private); `active_count` is a load-bearing swap-remove invariant an external write would desync. *Fix:* make the pool/counters private with accessors, exposing only the genuine tunables.

32. **[hardware] `hardware/pov_segmented.h:294-299 vs 321-323`** — `run_show()` mixes `hs::disable_interrupts()/enable_interrupts()` (platform wrappers) in the effect-publish bracket with raw `__disable_irq()/__enable_irq()` in the telemetry bracket 20 lines down. Correct on device, but bypasses the abstraction used immediately above. *Fix:* use the `hs::` wrappers in both brackets.

33. **[engine] `core/engine/transformers.h:320`** — `cos_threshold_min`/`cos_threshold_max` invert intuitive cosine ordering (`min` holds the *larger* cosine / nearest angle); the comparison is correct but reads backwards and a maintainer could "fix" it and silently break the ripple fast-reject band (no test would catch a perf-only regression). *Fix:* rename to `cos_near`/`cos_far`, or add a clarifying comment at the comparison site.

34. **[daydream] `daydream/tools/slider.js:35`** — `createSlider` returns `null` on a missing container while siblings (`initScene`, `formatFloatCpp`, codegen helpers) throw named errors; a caller that forgets the null-check gets an opaque `Cannot read properties of null`. *Fix:* throw a named `container #id not found` error for parity, or document the soft return.

35. **[animation] `core/animation/sprites.h:112`** — The pause-gate member is `paused_` in `Sprite` but `paused` in `Mutation`/`Driver`/`Lerp`. *Fix:* rename to `paused` for module consistency.

36. **[mesh] `core/mesh/mesh.h:561`** — Topology hashing quantizes interior angles with `std::round(ang*180/PI)` to whole degrees, so two congruent faces whose angles straddle a 0.5° rounding boundary (float noise) can land in different `topo_id` bins and never merge — smaller classes, reduced LUT coverage (never incorrect output). *Fix:* if class-merge robustness matters, compare on a coarser or hysteretic angle key.

_Cosmetic and documentation-only:_

37. ✅ **[math] `core/math/3dmath.h:726-727`** — `inv_stereo` comment claims `|z| >= 5000` "catches ... any point within ~1.1° of the pole," but `|z| = cot(ε/2)` gives ε ≈ 0.023°, not 1.1° (off by ~50×). *Fix:* change "~1.1°" to "~0.02°".

38. ✅ **[color] `core/color/color.h:1536-1544`** — `ProceduralPalette::get` clamps each channel to `[0,1]` then passes it to `srgb_to_linear_interp`, which re-clamps internally — a redundant per-channel clamp on the cold/bake path. *Fix:* drop the outer clamp, or comment the intent.

39. ✅ **[render] `core/render/scan.h:263-279`** — `BoundingSphere` ctor computes the full `vector_to_pixel(center)` for `center_theta` but discards `center_px.y` and recomputes `center_phi` from `center.y`, wasting a `phi_to_y` (cold path, negligible). *Fix:* use a theta-only projection helper.

40. ✅ **[daydream] `daydream/driver.js:511-539,895-903`** — `refreshLabels` allocates a fresh `labels` array and spreads `getLabels()` every rendered frame, and `coordsLabel` allocates a new `THREE.Vector3` per call, defeating the CSS2DObject pooling's intent. *Fix:* reuse a persisted labels array (`length = 0` + push) and write `coordsLabel` into a caller-supplied out-vector (mirroring `pixelToSpherical(out)`). *(Done: `refreshLabels` now reuses a persisted scratch array; `coordsLabel` has no in-repo caller, so its per-call allocation is off any render path and was left as-is.)*

41. ✅ **[daydream] `daydream/segment_controller.js:645-646`** — `renderParallel()` resets `timings`/`renderUs`/`arenas` each dispatch but not `results`, and `updateStats()` derives the per-segment "Range" cell from `results[s]`, so a fenced/not-yet-reported frame shows the prior generation's rectangle while the other columns show 0 (cosmetic). *Fix:* clear `results` alongside the other per-segment arrays at dispatch, or gate the range cell on `frameSeen[s]`. *(Done via the `frameSeen[s]` gate; clearing `results` at dispatch would break the overrun re-blit that re-uses the prior frame's `results`.)*

42. ✅ **[daydream] `daydream/sidebar_logic.js:89-92`** — `scrollArrowState` uses `deadzone = min(4, maxScroll/2)`, so at exactly `scrollLeft === maxScroll/2` both arrows hide even though content overflows both directions, and the two can never show simultaneously in the small-overflow regime (undercutting the header comment). *Fix:* use a strict deadzone only for the right edge, or `min(4, (maxScroll-1)/2)`.

43. ✅ **[tests] `scripts/wasm_smoke.mjs:259`** — Comment reads "MESHOP_1F(truncate)" but `truncate` is bound via the `_OP1U` [0,1]-clamped tier, not `_OP1F`. *Fix:* relabel to `MESHOP_1U(truncate)`.

---

## Closing Assessment

Nothing here calls for structural rework. The two Medium items are worth landing promptly (a real visual seam and a real accessibility defect); the Low items are a clean, well-scoped punch-list of consistency, guard-hardening, telemetry, and documentation polish that can be worked at leisure without risk. The eligibility bar was applied strictly — many plausible-looking candidates were examined and rejected as by-design (e.g. the antipodal slerp gate, the 2-gon Conway self-pair, composition polarity, ShapeShifter's negative radius, the AntiAlias non-compensation, the drop-on-overrun DMA policy), and those rejections are themselves a signal of how much prior care this codebase has absorbed.
