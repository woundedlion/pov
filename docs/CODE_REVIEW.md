# Holosphere / daydream — Code Quality Review

**Date:** 2026-06-29
**Reviewer:** Multi-agent audit (21 component reviewers + 4 independent adversarial validators)
**Commit reviewed:** `f161475b` (master)

## Scope

In scope: the C++ rendering engine (`core/`, except generated/third-party files), the
visual effects (`effects/`), the hardware drivers (`hardware/`), the WASM/Phantasm
targets (`targets/wasm/`, `targets/Phantasm/`), the native test suite (`tests/`), the
build/CI tooling (`CMakeLists.txt`, `scripts/`), and the daydream web simulator
(`c:\work\daydream`: app JS, segmented-POV pipeline, recorder, geometry tools).

Explicitly **out of scope** (per the review brief): `core/effects_legacy.h`,
`targets/Holosphere/Holosphere.ino`, `core/rotate.h`. Generated tables
(`core/color_luts.h`, `core/reaction_graph.cpp`) and third-party single-headers
(`core/FastNoiseLite*`) were reviewed only at their generator/integration boundary,
not as table data. Vendored `daydream/three.js/` and `node_modules/` were excluded.

## Methodology

Each component was audited end-to-end by a dedicated sub-agent that first read the
README to place its files in the whole architecture, then read every in-scope line and
cross-referenced real call sites before reporting. Every candidate finding was then
handed to a separate adversarial validator instructed to **refute** it — re-deriving
from the source, checking it was not correct-by-design, and confirming the proposed
fix. Only findings that survived independent validation appear below. Findings are held
to the `code-review-fix` eligibility bar: a real, demonstrable defect with a concrete,
minimal, non-performance-sacrificing fix.

---

## Executive Summary

This is an exceptionally mature, rigorously engineered codebase — among the best
hobby/art-engineering projects of its kind. The C++ engine core (math, color, memory,
rasterizers, mesh/Conway operators, animation, filters, hardware sync) is essentially
defect-free at the correctness and memory-safety level: eleven of the twelve
engine-core audits returned **zero** eligible findings after adversarial validation,
and the defenses are not accidental — invariants are pinned with `static_assert`,
guarded with always-on `HS_CHECK` traps, and documented at the point of use with the
*reason* a naive change would break them.

The defects that do exist are concentrated in the **periphery** — the JavaScript
simulator's segmented-render fault path, a single effect's long-uptime numerical drift,
documentation/comment drift, and test-coverage gaps for a handful of animation classes
and one CI self-check. There are **no Critical or High-severity defects**: the most
serious confirmed item (a controller re-dispatch after a defensive fault latch) is
Medium and self-recovers, and the rest are Low/Nit.

### Overall grades

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness** | A | Engine core is effectively flawless; one Medium JS-pipeline robustness bug and one Medium long-uptime numerical drift in a single effect. |
| **Memory safety** | A+ | Arena discipline, exhaustive bounds reasoning, fail-fast traps; no eligible defect found anywhere in scope. |
| **Concurrency / real-time** | A | Phantasm 1-wire sync protocol and ISR double-buffer are A+ and test-pinned; the lone gap is a JS worker-pool fault re-check (F-1). |
| **Error handling / robustness** | A | Principled fail-fast-vs-soft-degrade doctrine on device; a few Low JS/CI error-path gaps (F-9, F-10). |
| **API & interface expressiveness** | A | Crisp owned/borrowed and borrow/store type distinctions, compile-time domain coercion, single-source X-macro rosters; one dead API surface (F-4). |
| **Architectural elegance** | A+ | Compile-time resolution templating, variadic filter pipeline with `if constexpr` domain lifting, partitioned arena, and position-from-time flywheel are genuinely elegant. |
| **Performance / efficiency** | A+ | Hot/cold path discipline, LUT-backed math, zero-copy WASM views, branchless wrap guards; tradeoffs are measured and documented. |
| **Maintainability** | A | Comments are dense but load-bearing and accurate; the few defects here are stale comments/docs (F-7, F-12, F-13, F-14). |
| **Documentation** | A− | In-code doxygen is exceptional; most of this review's findings are README/comment drift, which is what docks it. |
| **Test coverage & rigor** | A− | Death harness, sim↔device parity recompiles, POV tiling proofs, and determinism scrambling are exemplary; gaps remain for several animation classes (F-3), one concept (F-15), and a CI self-check (F-11). |
| **Build / CI / tooling** | A− | Solid presets, install provenance, anti-drift roster checks; minor script-robustness and one ctest gap (F-9, F-10, F-11). |
| **Portability / determinism** | A+ | Bit-identical sim↔device determinism is treated as a first-class, tested invariant. |

### Per-component summary

| Component | Verdict |
|---|---|
| `core/` math, color, memory, concepts | A+ — zero eligible findings; impeccable edge-case discipline |
| `core/` sdf, scan, canvas, plot, filter, transformers | A+ — zero eligible findings (one stale README line, F-7) |
| `core/` animation, mesh, conway, solids, hankin, spatial, generators | A+ — zero eligible findings; manifold/Euler-verified |
| `core/` engine, platform, registry | A — zero eligible findings; anti-drift registry |
| `hardware/` sync, drivers | A — two documentation nits (F-12, F-13); load-bearing concurrency is correct & tested |
| `effects/` (27) | A — one Medium (F-2) + one Nit (F-14) across 27 effects |
| `targets/wasm`, `Phantasm` | A — one dead-API Low (F-4) |
| `tests/` (C++) | A− — coverage gaps (F-3, F-15) on an otherwise exemplary harness |
| daydream core JS | A — a cluster of Low/Nit findings (F-4, F-5, F-6, F-17) |
| daydream segmented + recorder | A− — one Medium (F-1) |
| daydream tools + build/CI | A− — F-9, F-10, F-11, F-16 |

---

## Prioritized Findings

No Critical or High-severity findings survived validation. Items are numbered
sequentially across all priority sections.

### Critical

*None.*

### High

*None.*

### Medium

1. **`SegmentController.tick()` re-dispatches a render after `composite()` latches a fault — doomed broadcast + leaked promise.**
   - File: `daydream/segment_controller.js:758-769` (root); interacts with `:534`, `:543` (composite faults), `:388-403` (`onWorkerFault`), `:228` (frame handler), `:749` (top guard)
   - Severity: Medium · Dimension: Concurrency · Confidence: High (control flow independently re-derived by validator)
   - Evidence: `composite()` (called only from `tick()`, `:759`) can call `onWorkerFault()` from its bounds/length pre-pass; `onWorkerFault` sets `faulted=true` and `renderInFlight=false`. Control returns to `tick()` with **no `faulted` re-check**, so `if (!this.renderInFlight)` is true and a fresh `render` is broadcast to the just-faulted pool. The frame handler early-returns while `faulted` (`:228`), so `pending` never decrements, `frameResolve` never fires (leaked promise), and `renderInFlight` stays stuck `true`. The faults are defensive invariant checks (reachable only if a fence/layout bug already let a bad rect through) and user-facing recovery still works via `create()→destroy()`, which is why this is Medium and not High. Untested: existing fault tests call `composite()` directly, never through `tick()`.
   - Fix: after the `pendingFrame`/`composite()` block in `tick()`, before the dispatch, add `if (this.faulted) return;`. No-op on the clean path. Add a test that delivers a fence-escaping out-of-bounds result through `tick()` and asserts no second `render` is broadcast.

2. **`Liquid2D` noise-time accumulator grows unbounded and eventually freezes the warp.**
   - File: `effects/Liquid2D.h:60-61, 107, 156-157, 274`
   - Severity: Medium · Dimension: Correctness · Confidence: Medium
   - Evidence: `accumulated_time` is driven by `Animation::Driver(..., wrap=false)` and passed directly as the FastNoiseLite Z-coordinate (`t * 0.5f`). With `Time Speed` up to 5.0/frame it climbs without bound; once its float ULP exceeds the per-frame increment the noise input is effectively constant and the warp visibly stalls/steps. The trig phases are already wrapped (`:117-118`); only the noise axis is exposed. The sibling effect `Flyby` guards the identical pattern (`Flyby.h:109`, `TIME_PERIOD = 65536.0f`), demonstrating the in-comment claim that the axis "cannot be wrapped" is not true for aperiodic OpenSimplex2.
   - Fix: wrap the accumulator as Flyby does — `accumulated_time = fmodf(accumulated_time, TIME_PERIOD)` at a large period (one seam per period is invisible for OpenSimplex2). Choose the period deliberately given Liquid2D's `t * 0.5f` halving doubles the effective noise-Z period.

3. **Several live `Animation::` lifecycle classes have no isolated unit test.**
   - File: `tests/test_animation.h` (gap); classes in `core/animation.h`: `ColorWipe` (`:1742`), `MobiusWarp` (`:1852`), `MobiusWarpCircular` (`:1906`), `MobiusWarpEvolving` (`:2095`), `Animation::Ripple` (`:2188`), `Animation::Noise` (`:2287`)
   - Severity: Medium · Dimension: Test coverage · Confidence: High
   - Evidence: `test_animation.h` directly tests the step/easing/`done()`/repeat lifecycle of essentially every other animation type, but these six are exercised only indirectly through the effects smoke pass (which cannot see an easing-application slip, snapshot-on-first-step regression, or wrong `done()` boundary). Note: `test_transformers.h` covers the pure `ripple_transform()`/`noise_transform()` *functions*, which are distinct types from the `Animation::Ripple`/`Animation::Noise` *classes* — the classes remain untested.
   - Fix: add lifecycle tests mirroring the existing `Transition` tests. `ColorWipe` is the highest-value/cheapest (assert the palette lerps from `from_snap` to `to_snap` and `done()` only on the final frame).

### Low

4. **`getParamGeneration()` is a dead contract — exposed and documented as a re-fetch protocol, but no consumer honors it.**
   - File: `targets/wasm/wasm.cpp:577-588` (getter + contract doc), `:1291` (embind export), `README.md:1872`; JS side `daydream/daydream.js:104-126, 226, 449-455`
   - Severity: Low · Dimension: API design · Confidence: High
   - Evidence: The bridge documents that "a consumer records this value … and re-fetches the definitions whenever it changes," but the only references in either repo are the C++ getter and the README — daydream never reads it. The safety the contract promises is actually delivered by construction: `applyEffect` rebuilds the param name list on every effect/resolution switch, and `syncGUI` binds by name with a `Math.min(names.length, values.length)` clamp. The counter `paramGeneration_` is still updated, so removing only the getter leaves a harmless write-only field.
   - Fix: prefer softening the doc to "provided for a *positional* (index-keyed) consumer; daydream binds by name and does not use this," or remove the getter + embind line + README row.

5. **`applyResolution` desyncs UI/URL from the engine if `setResolution` is rejected.**
   - File: `daydream/daydream.js:404-415` (handler), `:449-455` (subscriber wiring)
   - Severity: Low · Dimension: Error handling · Confidence: High (latent; unreachable today)
   - Evidence: `applyResolution` runs as the `appState` `resolution` subscriber, so `appState`/URL/dropdown already hold the new value by the time it runs; on `setResolution(...) === false` it logs and returns, leaving the engine on the old resolution while the UI/URL advertise the new one. Unreachable today because the two WASM-supported resolutions (`wasm.cpp:221-223`) exactly match the two daydream presets, so the call never rejects — but it is a latent footgun for any future build-gated third resolution (the documented extension point).
   - Fix: validate against `getSupportedResolutions()` (exists, `wasm.cpp:716`) before committing the `appState` write, or revert `appState` to the previous resolution on rejection.

6. **Deep-linked animated parameters are immediately overwritten by the animation (no auto-pause on hydration).**
   - File: `daydream/daydream.js:348-355` (hydration), `:357-367` (manual `onChange` auto-pause), `:104-126` (syncGUI snap-back)
   - Severity: Low · Dimension: Correctness · Confidence: Medium (intent-dependent)
   - Evidence: The manual edit path auto-pauses animation when an `animated` slider is touched; the deep-link hydration path pushes `?Param=value` straight through `setParameter` and never pauses, so for an animated param the engine re-drives the value next frame and `syncGUI` snaps the slider back. A shared/bookmarked animated-param value therefore does not "stick," asymmetric with a manual edit.
   - Fix: if animated deep links should hold, pause on hydrating any `p.animated` param; if the current behavior is intended, document it. (Owner intent call — not a unilateral fix.)

7. **README §7.2 describes the superseded `sqrt(1 - y²)`-only curve step model.**
   - File: `README.md:736` vs `core/plot.h:471-490` (`screen_step`), `:644-649` (comment)
   - Severity: Low · Dimension: Documentation · Confidence: High
   - Evidence: §7.2 states step size "is scaled by `sqrt(1 - y²)` — the sine of the polar angle." The live rasterizer instead sizes each sub-step from the full 2-D screen speed `sqrt(vx²+vy²)` (both `dlon_ds` and `dphi_ds`); the code comment explicitly contrasts with "the old sin(φ) longitudinal-only proxy." Per the repo's describe-current-state-only convention, the README is stale.
   - Fix: reword §7.2 to describe the screen-velocity sampler (each sub-step ~1 px apart in screen space, with a pole-oversampling floor). No code change.

8. **`Animation::RandomWalk` has no direct test of its on-sphere invariants.**
   - File: `tests/` (gap); `core/animation.h:1621`
   - Severity: Low · Dimension: Test coverage · Confidence: High
   - Evidence: `RandomWalk<W,CAP>` has deterministically-testable pure logic (smoothed angular-velocity integrator, on-sphere unit-length walk, degenerate-axis `normalized_or` fallback at `:1715`) but is exercised only inside the ~19 effects that use it. The unit-length invariant and the fallback branch are RNG-independent under a fixed seed.
   - Fix: add a seeded test that steps a `RandomWalk` ~50 frames and asserts the oriented vector stays unit-length each frame and total travel > 0; optionally force the degenerate case.

9. **`capture_screenshots.mjs` makes every effect report FAILED when `MAX_ATTEMPTS=0`.**
   - File: `scripts/capture_screenshots.mjs:22-27` (`numEnv`), `:41`, `:209-218`, `:227-230`
   - Severity: Low · Dimension: Error handling · Confidence: High
   - Evidence: `numEnv('MAX_ATTEMPTS', 6)` returns `0` for the input `"0"` (zero is finite and ≥ 0), so the grab loop never runs, `best` stays `null`, and `best.split(',', 2)[1]` throws. The TypeError *is* caught by the per-effect `try/catch` (`:227`), so it does not crash the run — but every effect is then logged FAILED, silently producing zero screenshots.
   - Fix: `const MAX_ATTEMPTS = Math.max(1, numEnv('MAX_ATTEMPTS', 6));` (or reject `< 1` in `numEnv`).

10. **`wasm_smoke.mjs` dereferences `m.stack` unguarded — an uncaught crash if the region is ever absent.**
    - File: `scripts/wasm_smoke.mjs:126-127` (and the same bare deref at `:198`)
    - Severity: Low · Dimension: Error handling · Confidence: Medium
    - Evidence: After iterating `Object.keys(m)`, the script does `const stack = m.stack; if (stack.high_water_mark === 0) …`. The enclosing block is `try { … } finally { engine.delete(); }` with **no `catch`**, so if `getArenaMetrics()` ever omits `stack`, the TypeError escapes the smoke harness instead of being counted via `fail()`. Defensive (the region is contractually present today), but cheap to harden, and the same pattern repeats at `:198`.
    - Fix: `if (!stack) fail(...)` before the dereference at both sites, or check the `stack` region presence once as a precondition.

11. **The LUT generator's own `--check` self-validation is not wired into the ctest suite.**
    - File: `tests/check_color_luts.cmake:12-19`, `tests/CMakeLists.txt:113-122`, `scripts/generate_luts.py:124-149`
    - Severity: Low · Dimension: Test coverage · Confidence: High
    - Evidence: `check_color_luts.cmake` runs the generator without `--check` and only token-diffs the committed table against fresh output — which would still pass if a libm/platform shift changed the table monotonically in *both* the committed file and the regen. The generator's `--check` path (monotonicity + sRGB round-trip ≤ 1 code) exists but runs only as a CI shell step, not via ctest, so the Windows ctest path skips it.
    - Fix: add a ctest mirroring `unit_color_luts` that runs `python generate_luts.py --check` (skip-coded when Python is absent).

### Nit / Polish

12. **README Phantasm pin table omits ID2 (pin 23) that the driver declares.**
    - File: `README.md:96` vs `hardware/pov_segmented.h:101-106` (`PIN_ID2 = 23`, `kIdStraps`), `:99` (doc)
    - Severity: Nit · Dimension: Documentation · Confidence: High
    - Evidence: The table says "ID: pins 21–22"; the driver declares three ID straps (21/22/23) and reads pin 23 only at `N = 8` (it is reserved-but-unread at the shipped `N = 4`). The README is correct for the shipped wiring but understates the declared hardware contract.
    - Fix: note that pin 23 is reserved for ID2 (used only at `N = 8`), or leave as-is matching shipped wiring.

13. **`hd107s_frame.h` comment claims the brightness byte is primed once, but it is re-written every pixel every frame.**
    - File: `hardware/hd107s_frame.h:97` (comment), `:173`, `:180`, `:216`
    - Severity: Nit · Dimension: Documentation / Maintainability · Confidence: High
    - Evidence: The class doc says the per-pixel `0xFF` brightness byte is primed at construction and "Only the B/G/R color bytes change at runtime," but `packPixel()` and `load()` unconditionally re-store `dest[0] = 0xFF` on every pixel, making the constructor priming dead for image frames.
    - Fix: correct the comment (low-risk), or drop the redundant per-pixel store as a micro-optimization — but the latter must be validated against the tests since `load()` only touches image slots, never the trailing black frame.

14. **`MindSplatter` memory comment cites the wrong device arena size.**
    - File: `effects/MindSplatter.h:124-125` vs `core/memory.h:37` (`DEVICE_GLOBAL_ARENA_SIZE = 330 * 1024`)
    - Severity: Nit · Dimension: Documentation · Confidence: High
    - Evidence: The comment says "316 KiB of the 335 KiB device arena, leaving ~8 KiB"; the arena is 330 KiB, not 335. The binding `static_assert` uses the correct literal, so the build is unaffected — comment-only. (The "~8 KiB" is gross slack after the scratch carve; ~2 KiB remains after the 6 KiB aux reserve.)
    - Fix: change "335 KiB" → "330 KiB" and state both the gross (~8 KiB after scratch) and net (~2 KiB after aux) figures to avoid ambiguity.

15. **The `Tweenable` concept has no compile-time pinning test.**
    - File: `tests/test_concepts.h` (gap); `core/concepts.h:323-327`
    - Severity: Nit · Dimension: Test coverage · Confidence: High
    - Evidence: `test_concepts.h` thoroughly covers `FunctionRef`/`Fn`, and `StoredFunctionRef`'s rvalue-rejection is already pinned (`test_concepts.h:118-119`) and `PipelineRef` is covered (`test_canvas.h:571`) — but the `Tweenable` concept that constrains `Lerp`/`Transition` subjects has no `static_assert` confirming a conforming type satisfies it and a non-conforming one does not.
    - Fix: add one `static_assert` pair (e.g. `Tweenable<Lerpable>` holds, `Tweenable<int>` does not). Compile-time only.

16. **`palette_math.js` `SPLIT_COMPLEMENTARY` case declares `const` without a block (`no-case-declarations`).**
    - File: `daydream/tools/palette_math.js:309-313` (vs the correctly-braced `ANALOGOUS` case at `:318-323`)
    - Severity: Nit · Dimension: Style/Convention · Confidence: High
    - Evidence: `case "SPLIT_COMPLEMENTARY":` is immediately followed by `const complement = …` with no braces; the identifier leaks into the shared switch scope. Behavior is correct today but it trips the `no-case-declarations` lint and is a reorder hazard; the sibling case already uses braces.
    - Fix: wrap the case body in `{ }` to match `ANALOGOUS`.

17. **In segmented mode the main-engine arena metrics are recomputed and written into hidden DOM cells every frame.**
    - File: `daydream/daydream.js:512-513`, `daydream/driver.js:791-806`
    - Severity: Nit · Dimension: Performance · Confidence: High
    - Evidence: `getArenaMetrics()` always returns the (idle) main-thread engine's metrics, and `driver.updateStats` writes them into the global stats cells each frame. In segmented mode those cells are hidden (`#global-stats-desktop`/`#stats-bar` are display:none; the visible HUD is the correct per-segment `#segment-stats` table sourced from the workers). So there is **no stale-display bug** — only a small redundant `getArenaMetrics()` call + text writes into invisible nodes.
    - Fix: skip the arena block in `driver.updateStats` when `segments.active` (efficiency only).

---

## Investigated and rejected (not findings)

For transparency, two candidates were raised by component reviewers and then **refuted**
by independent validation:

- **`solid_codegen.js` snub `twist` suffix collision (doc gap)** — *Rejected.* The cited
  example is arithmetically wrong (`pctSuffix(0.144)` → `"14"`, `pctSuffix(0.145)` →
  `"15"`; no collision), and the `pctSuffix` docblock (`tools/solid_codegen.js:56-62`)
  already documents the rounding-granularity collision risk for *all* fractional params,
  not just `t`. The collision class is real and already disclosed; no fix needed.
- **Segmented-mode arena HUD shows "stale/unrelated" figures** — *Reframed, not a display
  bug.* The visible segmented-mode HUD is the per-segment table, which is correct; the
  main-engine write lands in hidden cells. The residual (a redundant per-frame write) is
  captured as the Nit F-17 above.

---

## What is genuinely excellent

- **Fail-fast as an enforced doctrine, not a slogan.** Always-on `HS_CHECK` traps guard
  cold seams (arena OOM, capacity/bounds, non-manifold input) while per-pixel hot paths
  use stripped `assert` backed by a cold bind-site trap — and a death harness asserts the
  exact `SIGILL`/`STATUS_ILLEGAL_INSTRUCTION` shape so the safety net is *tested*, not
  assumed.
- **Compile-time invariants over runtime hope.** Span-count CSG algebra, effect-roster
  and resolution X-macros, param-marshal ordering, preset field counts, and index-width
  ceilings are all `static_assert`-pinned, converting whole classes of drift into build
  errors.
- **The Phantasm position-from-time flywheel.** Deriving column position from the
  free-running cycle counter (`x = f(now − epoch)`) makes masked-IRQ windows quantize lag
  but structurally *incapable* of dropping columns — and the whole 1-wire protocol is
  host-testable with a multi-board drift simulator.
- **Sim↔device determinism as a first-class, tested invariant.** Shared PRNG seeds,
  uint32-narrowed time, and extra test TUs that recompile the engine under the *shipping*
  fast-math flags and device `H_OFFSET` keep the WASM simulator bit-faithful to hardware.
- **Memory architecture.** The partitioned arena (persistent + two scratch pools) with
  RAII `ScratchScope`/`Persist`, dual-stamp dangling detection, and owned-vs-borrowed
  type distinctions is a genuinely sophisticated solution to deterministic memory on a
  fragmenting embedded heap.
- **Documentation that pre-litigates correctness.** In-code comments routinely explain
  *why* a non-obvious guard exists and what a naive "cleanup" would break — which is what
  made an adversarial review of this depth tractable.
