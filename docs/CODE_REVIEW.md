# Holosphere / Daydream — Comprehensive Code Quality Review

**Date:** 2026-06-22 (fourteenth audit; supersedes the thirteenth. Its fix list 1–27 is essentially fully resolved in the tree; the two still-open carry-overs and the one standing rejection are reconciled in §5.)
**Tree:** Holosphere `master`; daydream reviewed at its current `master`.
**Scope:** the entire two-repo product — the Holosphere C++ engine + firmware (`core/`, `effects/`, `hardware/`, `targets/wasm`, `tests/`, `scripts/`, `tools/`, build) and the daydream web simulator (`*.js`, `tools/*.js`, `tests/*.test.js`). Out of scope per instruction: `core/effects_legacy.h`, the `*.ino` target sketches, and `core/rotate.h`. `core/FastNoiseLite.h` (vendored) and the generated data tables (`core/color_luts.h`, `core/reaction_graph.cpp`) were assessed for integration/provenance only.
**Method:** the README architecture was read in full, then ~18 independent sub-agents each audited a coherent component, reading every in-scope file with the README as shared architectural context and grading objective quality (correctness, UB, ISR/memory-ordering, bounds, overflow, performance, testability) and subjective quality (architectural elegance, interface expressiveness, idiom). The prior `CODE_REVIEW.md` was **not** consulted during examination; it was reconciled only afterward (§5).

---

## 1. Executive Summary

Holosphere is a persistence-of-vision LED sphere rendering engine and its bit-identical
WebAssembly simulator. The same header-only C++17 engine — every rendering class templated
on `<int W, int H>` — compiles to Teensy 4.0 firmware (single-board Holosphere and
4-board synchronized Phantasm) and, via Emscripten, to the browser simulator wrapped in
Three.js.

Read cold against the current tree by eighteen reviewers spanning both repositories, the
product holds the high bar its audit series has established. **Across well over a hundred
source files, not a single Critical or High-severity correctness, memory-safety, or
concurrency defect was substantiated.** The recurring finding was the *absence* of the
defects one expects: the zero-copy WASM view-detachment hazard is defended at every
consumption site; the multi-board sync protocol's failure mode is provably "missed, never
wrong"; arena lifetimes are explicit and trap-guarded; numeric casts are uniformly clamped
against NaN/UB; and the fail-fast (`HS_CHECK`) discipline is applied exactly where the design
philosophy prescribes (cold paths) and deliberately withheld from hot paths. Every reviewer
of the C++ core independently reported the same pattern: the residual defects are
overwhelmingly **internal inconsistencies with the codebase's own exacting standards** —
the one `acosf` argument that isn't clamped where every sibling clamps, the one fade multiply
that truncates where every sibling rounds, the one comment copied from a sibling effect that
no longer matches its file. When a codebase's flaws are best described as exceptions to its
own rules, the rules have largely won.

The findings that remain concentrate in four bands, none load-bearing for correctness today:

1. **Test-coverage breadth on the JS browser-glue layer** — `driver.js` (~880 lines),
   `daydream.js` (~805 lines), `sidebar.js`, and `recorder.js`'s session lifecycle are
   essentially untested — and **behavioral blindness of the C++ effect smoke harness**
   (~22/28 effects are only smoke + determinism checked, so a wrong-but-non-black,
   reproducible regression passes).
2. **Device-budget items that cannot be verified in CI** because the Teensy build compiles
   nowhere automated (e.g. `SplineFlow`'s 240 KB trail buffer on a 335 KB arena).
3. **A small number of latent edge cases worth a confirming test** — chiefly in the Phantasm
   sync timing math (the beacon mod-64 resync precondition and the same-tick burst/fold
   ordering) and a couple of unclamped transcendental arguments.
4. **Comment/naming/consistency polish**, including two violations of the project's own
   "no finding-number references in comments" convention.

**The grade is A−.** As in every prior audit, the distance to a flat A lives almost entirely
at the **test-coverage frontier** — not the test *infrastructure*, which is A+ (a forked
death harness proving every `HS_CHECK` trap fires with the exact illegal-instruction signal;
a fault-injecting four-board sync simulator; injected-clock cross-run determinism; analytic
oracles over Euler characteristic, KD-tree brute force, and SDF cull-conservativeness) — but
the test *coverage* of the most delicate paths: the untested JS render/orchestration glue and
the behaviorally-unpinned effect roster. The engine, color, memory, mesh, hardware-protocol,
WASM-bridge, build/CI, and test-infrastructure subsystems each grade A or above on their own
terms.

---

## 2. Quality Dimensions — Letter Grades

| Dimension | Grade | Rationale |
|---|---|---|
| **Architectural elegance** | **A+** | Compile-time `<W,H>` specialization, the World/Screen/Pixel filter pipeline with automatic compile-time domain conversion, the state-mutation-vs-render separation in the animation system, the SDF/scan split, the CRTP reaction-diffusion base, and the three-layer flywheel sync are genuine design contributions. Subsystems compose rather than entangle. |
| **Interface expressiveness / API design** | **A** | `StaticPalette` compile-time modifier chains, `SolidBuilder` fluent Conway chaining, the `generate()` arena wrapper, `FunctionRef`/`inplace_function` zero-alloc callables, deleted rvalue overloads enforcing borrow contracts, and `registerParam()` reflection are expressive and hard to misuse. Minor sharp edges (`ParamDef::set` raw-write bypass; public mutable `Presets` indices). |
| **Correctness & robustness** | **A−** | No substantiated Critical/High bugs across 100+ files; the degenerate-geometry policy (`normalized()` trap vs `normalized_or()` soft-degrade) is principled and applied per call site. Held below A by a scattering of latent edges (unclamped `acosf` in `fib_spiral`, the two sync-timing items, the Volume halo early-out) whose siblings guard. |
| **Memory safety** | **A** | Arena bounds math is wrap-safe (subtractive no-wrap checks + independent multiply-overflow trap) and alignment-correct; dual-stamp use-after-free detection on `ArenaVector`/`ArenaSpan`; deleted rvalue ctors prevent dangling borrows; zero per-frame heap on any path. |
| **Concurrency (ISR / sync / workers)** | **A** | Single-writer discipline, correct relaxed-atomic reasoning for single-core M7, position-from-time flywheel, airtight Transferable ownership and generation-fencing in the worker pool. |
| **Numerical robustness** | **A−** | Pole/antipode/normalization edges handled deliberately; every `acos`/division domain-clamped *except* one (`fib_spiral`); fast-math approximations carry measured error budgets and exact-trig escapes where needed. |
| **Error handling / fail-fast** | **A** | A faithful application of the doctrine: cold-path invariant violations trap via stdio-free `HS_CHECK` (survives `NDEBUG`, death-harness-verified); transient conditions get bounded soft handling; the WASM/JS boundary logs-and-rejects untrusted input rather than trapping. |
| **Performance** | **A−** | LUTs, interval culling, branchless ISRs, packed-SIMD blits, hot/cold path partitioning, zero per-frame heap. A few effects (`Dynamo` per-pixel band walk, `FlowField`/`Voronoi`/`Raymarch`/`ShapeShifter`) carry per-frame cost cliffs gated on sliders rather than the pixel budget — flagged honestly in-comment. |
| **Readability & naming** | **A−** | Consistently descriptive; comments are load-bearing, not decorative — the *why* sits at the line that depends on it. A few terse names (`G`, `ttls_`), over-long comment essays, and a misleading `size` field. |
| **Comments & documentation** | **A** | The 2,100-line README is a reference-class architecture document, accurate to code (palette lists, register conventions, sync constants verify). Inline Doxygen is near-universal and high-value; a small number of comments drift (noted per-finding). |
| **Testability** | **A−** | C++: outstanding (oracles, death harness, determinism pass, tiling proofs, parity TUs). JS: pure modules well-covered, browser glue not. |
| **Test coverage** | **B+** | The grade-limiter, as in every prior audit. C++ suite is **A** (~2,600 oracle-driven assertions); dragged down by untested JS browser-glue modules and behaviorally-blind effect smoke runs. |
| **Build / CI / tooling** | **A−** | Three-layer test gate, runtime WASM smoke, install-provenance trio, self-gating Teensy size/layout gate. Parser fragility against toolchain-output drift is the main risk. |
| **Consistency** | **A** | Idioms (arena params, `HS_CHECK`/`assert` split, `Filter::` namespacing, narrowing casts, X-macro rosters, Doxygen) are uniform across a large surface. |
| **Portability (Teensy / WASM / host)** | **A** | `platform.h` cleanly abstracts the targets with a *documented divergence ledger*; host-testable cores split out of Arduino-only shells; the NaN-clamp contract is guarded by a compile-time `#error`; deterministic seeded RNG gives bit-identical multi-target renders. |
| **OVERALL** | **A−** | Code quality is A/A+ across the engine; the composite is held at A− by the test-coverage frontier (untested JS render glue + ~22/28 behaviorally-unpinned effects) and a thin band of latent hazards. |

---

## 3. Notable Strengths

- **Provably-safe multi-board synchronization.** `pov_sync.h` derives column position from a
  free-running cycle counter (not interrupt counting), making masked-IRQ windows
  correctness-neutral by construction, paired with an odd-only distance-2 symbol codec whose
  worst case is a *missed* symbol (self-healed next half-rev), never a *misclassified* one.
- **Fail-fast that is itself tested.** `HS_CHECK` traps are verified by a cross-process death
  harness that asserts the exact `SIGILL`/`STATUS_ILLEGAL_INSTRUCTION` status and hard-fails
  CI if unrunnable — the safety net is proven, not assumed.
- **The hardest WASM hazard, handled everywhere.** The zero-copy `getPixels()` view detaches
  on heap growth; every consumption path re-fetches and the three-way buffer alias
  (`wasmMemoryView` ≡ `Daydream.pixels` ≡ `instanceColor.array`) is asserted on both the
  single-engine and segmented composite paths.
- **Numeric discipline.** 16-bit linear-light pipeline end-to-end, OKLCH perceptual palette
  interpolation with chroma-preserving gamut clipping, round-to-nearest fixed-point, and
  UB-safe float→int casts (every cast preceded by a NaN-folding clamp).
- **Memory model as a design feature.** A single 335 KB arena, repartitionable per-effect,
  with explicit `Arena&` parameters at every call site, RAII scratch scopes, `Persist<T>`
  evacuation, and dual-generation use-after-free detection.
- **Parity-by-delegation in the tools.** `splines.html`/`solids.html` and the generative
  palette tool route their math through the *same WASM engine code* rather than reimplementing
  it, eliminating the most likely source of sim-vs-device drift.
- **A reference-class README** that documents not just the *what* but the *why* of every major
  decision, down to AC timing tables for the 1-wire sync protocol.

---

## 4. Prioritized Fix List

Every actionable defect found this pass, grouped by severity and numbered sequentially (a
single global sequence, matching the codebase's finding-number convention). Each entry:
`file:line — problem → fix`. Items independently re-confirmed open from the thirteenth audit
are marked **(carry-over)**. A leading ✅ marks an item resolved during/after this review;
❌ marks a deliberately-rejected change preserved so it is not re-raised.

### Critical

*No Critical findings were substantiated. No correctness, memory-safety, concurrency, or
data-loss defect was found anywhere in scope.*

### High

1. `daydream/driver.js` (~880 lines) — the core render/Three.js/loop driver has **zero test coverage**: nothing pins its memory-view refresh, the buffer-alias contract `segment_controller` depends on, or frame cadence. → Extract the DOM-free logic and unit-test it with a stubbed Three.js.
2. `daydream/daydream.js` (~805 lines) — the main orchestrator (effect/param plumbing, resolution presets, state wiring) is entirely untested. → Add coverage for its pure orchestration paths behind driver/gui mocks.
3. `daydream/recorder.js` — only pure helpers are tested; the entire `MediaRecorder` session lifecycle (start/stop/captureFrame and the explicitly-reasoned rapid stop→start stale-session isolation) is untested — the module's subtlest, highest-risk logic. → Add a fake-MediaRecorder test driving start→capture→stop and a stop→start→old-onstop race.
4. `tests/test_effects.h` — **(carry-over, prior #10)** ~22 of 28 effects are covered only by `smoke_one` (no-crash) + `determinism_one` (reproducible); both are behavior-blind, so a wrong-but-non-black, reproducible regression passes. Only 6 effects have behavioral white-box oracles. → Add one pinned numeric property each for the highest-risk untested effects (Raymarch, HopfFibration, MobiusGrid, SphericalHarmonics, Voronoi).
5. `effects/SplineFlow.h:21,40` — `MAX_TRAILS = 30000` makes `Filter::World::Trails<W,30000>` allocate ~240 KB of the 335 KB device persistent arena, with no budget justification or `static_assert` (unlike `Dynamo`'s documented 80 KB cap), and the device build compiles nowhere in CI. → Add a `static_assert` against the device arena and a budget comment; confirm fit on real hardware.
6. ✅ `hardware/pov_sync.h:1445,1534` — the beacon mod-64 rev resync is correct only while the true slip is < 32 revs, an undocumented and unenforced dependency between `beacon_period_revs`, the join grid, and the 6-bit channel width. → Add a `Config::valid()` clause asserting `beacon_period_revs < 32` (or the actual derived bound) so the resync precondition is enforced rather than assumed.
7. ✅ `hardware/pov_sync.h:841-844,1145,1170` — `on_epoch_symbol`'s `j`-inference relies on the local boundary crossing having already folded `rev_in_effect`, but `tick()` handles the burst (line 1145) *before* the fold loop (line 1170); a burst arriving in the same tick as its boundary fold could infer `j` one short and commit one revolution late. → Add a host test for the same-tick burst-and-fold case; fold before handling the burst if it fails. *(Validated: no reorder needed — `handle_burst`'s backstop `apply_flip(ZERO)` already folds `rev_in_effect` before the `on_epoch_symbol` j-inference, so the same-tick case commits in lockstep; the later fold-loop flip is deduped. Added `test_epoch_same_tick_burst_fold` pinning the lockstep across every copy j, plus a comment marking the ordering as load-bearing.)*

### Medium

8. ✅ `core/platform.h:771,789,818` — `beatsin8`/`beatsin16`/`map8` compute `highest - lowest` as an unsigned subtraction with no guard, so a caller passing `lowest > highest` underflows and escapes the `[lowest, highest]` range their docstrings promise. → Document the precondition or clamp. *(Documented the `lowest <= highest` precondition on all three host mocks. Not clamped: real device `<FastLED.h>` subtracts identically and is unguarded, so clamping only on the host would break the sim/device parity contract.)*
9. ✅ `core/platform.h:1305,1347` — profiling accumulates `uint32_t` cycle deltas with no handling for the 600 MHz `CYCCNT` ~7 s wrap; a scope spanning a wrap underflows to a huge delta and `pct` can exceed 100 (cold path). → Accumulate in `uint64` or document that a scope must not span a counter wrap. *(`CycleCounter::cycles` and the `log_node` math widened to `uint64`, fixing the accumulation overflow that a multi-frame run hits in ~7 s. The per-scope delta stays a 32-bit subtraction matching the hardware CYCCNT register, with a documented precondition that no single scope spans a full wrap.)*
10. ✅ `core/platform.h:599` — `SerialMock::printf` uses `vsnprintf` (full float formatting) while the device path uses integer-only `vsniprintf`; a `%f` works on host but silently drops on device. → Flag the divergence as is already done for `check_fail`. *(Added a `@warning` documenting the divergence and that `%f`/`%g` must be avoided in device-bound messages; also corrected the stale "mirroring hs::log's vsnprintf path" line — hs::log uses vsniprintf. Not unified on vsniprintf: it is a newlib-only extension absent from the host toolchains, which is why check_fail already branches on `#ifdef ARDUINO`.)*
11. ✅ `core/concepts.h:22` — stale comment claims `std::function` is the WASM callable backend; the WASM path now uses `hs::inplace_function`. → Update.
12. ✅ `core/3dmath.h:638-642` — `Complex::operator/` uses ad-hoc magic thresholds `1e-6f`/`1e-12f` in a numerically-sensitive divide, contradicting the file's own "use the named `math::EPS_*` constants" directive. → Promote to named constants.
13. ✅ `core/geometry.h:518-523` — `fib_spiral` calls `acosf(1 - 2*(i+eps)/n)` with no clamp; at `i == n-1` floating error can push the argument below −1, making `acosf` return NaN. This is the one place an unclamped `acos` can NaN in a codebase that otherwise clamps every such call. → Clamp the argument to `[-1,1]`.
14. ✅ `core/color.h:132-140` — `Pixel16::operator*(float)` truncates (no `+0.5f`) while the rest of the file rounds-to-nearest; this is the per-frame fade multiply, so every faded channel is systematically biased down each frame — the exact dimming the file documents avoiding elsewhere. → Add the round or document why truncation is intended here.
15. ✅ `core/color.h:446-450` — `blend_alpha` truncates the lerp weight with a bare "intentional" comment giving no reason, diverging from `Color4::lerp`/`Gradient::get` which round the identical weight. → State the rationale or round for consistency.
16. ✅ `core/canvas.h:471,488` — `registerParam(float*)` captures the default as `*ptr` without clamping into `[min,max]`, but every later `updateParameter` clamps; an effect initialized out-of-range advertises an out-of-range default that snaps on first edit. → Clamp the captured default or add a cold-path `HS_CHECK(min <= *ptr <= max)`. *(Added the cold-path `HS_CHECK(*ptr in [min,max])` at the registration seam, alongside the existing capacity/duplicate/min<=max guards; NaN defaults trap too. The bool overload captures `(float)*ptr` ∈ {0,1} ⊂ [0,1], so it needs no guard.)*
17. ✅ `core/canvas.h:272-282,332-338` — `ParamDef::set()` plus a public mutable `find()` returning a writable `ParamDef*` lets callers bypass the readonly/finite/clamp gate in `updateParameter`; safety rests on caller discipline. → Return a const handle from the public lookup and confine mutation to a typed API. *(The untrusted WASM boundary was already const-gated — `getParameters()` returns `const ParamList&` and `updateParameter` is the sole write path. Closed the residual in-engine bypass: the mutable `find()`/`begin()`/`end()` are now `private` and reachable only by a `friend class Effect;`, so the public lookup hands every other caller a `const ParamDef*` that cannot call the non-const `set()`. No external/subclass caller used the mutable overloads.)*
18. ✅ `core/memory.h:78,470` — `Arena::allocate(0, align)` returns a non-null bump pointer without reserving storage, with no documented size==0 semantics; a future caller could treat the result as ownable. → Document or assert `size > 0`. *(Added `HS_CHECK(size > 0)` on the cold allocation path and documented the zero-size aliasing hazard; every real caller passes a positive byte count, and `ArenaVector::bind` already guards `exact_capacity > 0`, so no live path regresses.)*
19. ✅ `core/sdf.h:777-783,873-890` — `Union`/`SmoothUnion` interval generation is safe only by the `2×cap` `MergedIntervalBuffer`/`scan_region` sizing and overlap-collapse, but this reasoning is undocumented at the site and the `push_interval` overflow trap can fire on a row destined for full-scan fallback. → Document the buffer-sizing invariant; consider emitting `[0,W)` explicitly instead of returning false. *(Documented the invariant at both `get_horizontal_intervals` sites: each child caps at `kIntervalSpanCap` and `MergedIntervalBuffer` is `2 * kIntervalSpanCap`, so two full children fit and the trap cannot fire even though both push before the fallback check; the k-pad widens spans without adding any. Kept `return false` — the full-width scan it requests is equivalent to and cheaper than emitting `[0,W)` — and noted that.)*
20. ✅ `core/scan.h:1197-1210` — `Volume`'s past-the-back-of-sphere early-out breaks without updating `closest_d`, so the cull can never tighten an already-recorded AA halo and the occlusion probe can march from a point behind the sphere. Documented as a cosmetic limit — → confirm intended, or re-evaluate distance at the exit point before breaking. *(Confirmed intended. The exit point is past the back of the sphere, so its signed distance is large/positive and can never lower the running `closest_d` — re-evaluating there only wastes a `shape.distance()` call. And the occlusion probe already seeds from `closest_local` (the closest approach), not the terminal point, so the unrecorded exit sample doesn't weaken the halo cull. Added a note at the break site stating both; no code change.)*
21. ✅ `core/conway.h:1047-1052` — `relax` convergence is tested on the raw pre-normalize spring force, not the post-normalize angular move, so the comment's "largest per-vertex move below ~1e-4 rad" overstates precision (conservative, but inaccurate). → Reword or measure the post-normalize delta. *(Reworded the comment to state it tests the raw pre-normalize force magnitude, which is a tight conservative proxy for the post-normalize geodesic step (for a small near-tangential force the angular move is ~|force| rad). Left the convergence math unchanged — measuring the post-normalize delta would cost an extra pass for no behavioral gain.)*
22. ✅ `effects/HankinSolids.h:33,274-279` — the scratch-B 32 KB sizing comment accounts for the compaction phase, but `compile_hankin` also uses scratch-B as generation workspace; the 32 KB must cover the *larger* of the two phases. → Verify the peak high-water at H=144 and tighten the comment to name the dominating phase (device-only path). *(Tightened the comment to name both non-overlapping scratch-B peaks — generation (the solid's `.generate(a,b)` Conway intermediates + `classify_faces_by_topology` pairing scratch; note it is `generate_base_solid`, not `compile_hankin`, that uses `b`) and compaction (the morph Persist set CompiledHankin ~10 KB + MeshPaletteBank ~15 KB ≈ 25 KB) — and to flag that the H=144 high-water cannot be measured in CI (no automated Teensy build) and must be confirmed on hardware. Comment-only; the budget is unchanged.)*
23. ✅ `effects/SphericalHarmonics.h:305-307` — `transition = 0.03f` makes the "quintic-smoothed" crossfade a near-binary step for almost all field values; the comment oversells the smoothing. → Widen `transition` or correct the comment. *(Corrected the comment: it is a narrow quintic-smoothed anti-aliasing seam at the zero-crossing, not a wide crossfade — `transition` is the half-width of the AA band, so only |val| < 0.03 is on the ramp and the lobes meet at a crisp boundary by design (widening would blur the lobe edge). Kept 0.03 — the crisp boundary is the intended look, so this is a comment fix, not a visual change.)*
24. `effects/GSReactionDiffusion.h:135` — `STEPS_PER_FRAME = 16` (vs BZ's 2) is the single largest per-frame cost and reads as arbitrary. → Add a one-line justification for the 16× substep count.
25. `effects/FlowField.h:158` — the 600-particle × 3-noise-call per-frame force loop is the field group's second-heaviest cost and, unlike `Raymarch`, carries no device budget note. → Add one or confirm it fits the per-column budget.
26. `effects/Dynamo.h:189-252,305` — `color()` linearly walks every active palette boundary and is invoked per replayed trail point, giving `O(points × boundaries)` `angle_between`+LUT scans on the hot path. → Hoist the per-flush band lookup or precompute; profile on device.
27. `effects/MindSplatter.h:107-115` — the pool-footprint comment's hand arithmetic ("~19 KB headroom") is thin and rounds KiB/KB loosely while that margin must also cover scratch, palette, and emitters. → Verify actual `sizeof` and add a `static_assert` on the pool against the carved arena.
28. ✅ `targets/wasm/wasm.cpp:588-645` — `getParameterDefinitions()` and `getParamValues()` are each internally index-aligned but are two separate JS calls; the parameter *set* is frozen at `init()` today but not contractually, so a cached definitions array could mis-describe a later value stream. → Document the re-fetch contract or expose a monotonic param-generation counter. *(Added a `paramGeneration_` counter bumped on every set-changing transition — `setEffect()` install and `setResolution()` teardown, but not per-value `setParameter()` — exposed via a new `getParamGeneration()` binding. Documented the re-fetch contract on the accessor and cross-referenced it from both `getParameterDefinitions()` and `getParamValues()`: cache the counter with the descriptors and re-fetch when it changes.)*
29. ✅ `targets/wasm/wasm.cpp:308,313` — comment says the default effect is `IslamicStars` while the code installs `DistortedRing`. → Reconcile. *(Rewrote the comment to name both defaults and state they differ by design: `DistortedRing` is the C++ bootstrap default — any real registered name keeps the engine valid for the first instant/headless use — while daydream's JS frontend default `IslamicStars` (or `?effect=`) overrides it almost immediately. Verified daydream.js:177 does default to `IslamicStars`.)*
30. ✅ `tools/teensy_gate.py:93-96` — `_TS_REGION_RE` requires `\s{2,}free for`; the size-authoritative `teensy_size` output has variable spacing, so a future single-space variant silently yields "region-missing." → Loosen to `\s+free for` and add a real-capture regression test asserting all three regions parse. *(Loosened to `\s+`; documented why the lazy `.*?` + unique `free for` literal makes single-space safe without over-consuming the blob's internal spaces. Added `test_parse_teensy_size_single_space_before_free` asserting all three regions (FLASH/RAM1/RAM2) parse used+free from a single-space variant; the 2-space and real-capture fixtures remain covered by the existing parser/RealCapture tests.)*
31. ✅ `tools/teensy_gate.py:249-281` — the symbol-magnitude check iterates every name match including size-0 UND/duplicate rows, which could emit a spurious "symbol too small." → Filter to defined symbols (not UND / size>0) before the magnitude/region asserts. *(Filtered `matches` to `ndx != "UND" and size > 0`, so a same-named UND reference row no longer mis-derives a region from its null address or trips the magnitude floor; a name present only as UND still hard-fails `symbol-not-found`. Added `test_und_reference_row_does_not_trip_magnitude_or_region` — verified the injected UND row would otherwise fire both `symbol-too-small` and `symbol-wrong-region`.)*
32. ✅ `tools/teensy_gate_extra.py:39-43` — a toolchain step failing surfaces as an opaque traceback or an empty-dict "region-missing," misattributing a parser/toolchain break to a budget break. → Wrap in try/except emitting a distinct `::error::` annotation. *(Wrapped the ARM-tool capture/parse in a try/except for `OSError`/`subprocess.SubprocessError`, emitting a distinct `::error::teensy-gate: a toolchain step failed…` annotation and `exit(2)` (vs the `exit(1)` budget-violation path). Added an empty-regions guard emitting a separate `::error::` when no FLASH/RAM1/RAM2 parse out, so a format break is never reported as a budget failure.)*
33. `daydream.js:452-469` — `applyResolution()` can double-apply the effect: the forced `appState.set('effect', …)` fires `applyEffect()` synchronously *before* `updateResolution()` rebuilds the dot mesh/sidebar, so the effect GUI is constructed mid-resize. → Reorder so the effect switch happens after `updateResolution()`.
34. `daydream.js:457` — the resolution-preset field named `size` is actually the dot *radius*, inviting misreading as pixel count. → Rename to `dotSize`.
35. `tools/lissajous_math.js:40-46` — `lissajous()` omits the `.normalized()` the engine applies (`geometry.h:760`), so the JS preview diverges slightly from the device off the `sin(m2·t)≈0` points. → Apply the same normalization (a parity hazard for an export tool).
36. `tools/lissajous_math.js:54-83` — `findBestRationalRatio` only scans positive `M/N`, so a negative target ratio snaps to the wrong (sign-flipped) figure. → Take `abs` and reapply `sign`, or document non-negative inputs.
37. `daydream/tools/color.js` — `srgbToOklch`/`oklchToLinearRgb` are asserted only for finiteness/range, not against golden reference values, so a wrong-but-finite matrix or hue error passes. → Pin a golden OKLCH triple per function.
38. `daydream/tests/color_parity_wasm.test.js:89` — the palette parity tolerance is ±2 16-bit LUT steps, looser than the "within one step" the docstring claims; a 2-step systematic bias passes as parity. → Tighten to ±1 (or justify) and add a fixed golden value to catch a uniform offset.
39. `platformio.ini:69` vs `CMakeLists.txt:71` — the device firmware builds without `-ffast-math` while WASM builds with it, so simulator float results are not bit-identical to the device (presumably intended given the integer-wrap-guard parity contract). → Confirm acceptable and document where a reader of either build file would see it.
40. ✅ `tests/test_effects.h:175` — the smoke assertion only proves `get_pixel` is a stable pure accessor; it does not assert the frame is non-trivial, so a silent regression-to-black passes for effects expected to light up. → Track a "rendered non-black at least once across the sweep" assertion. *(Added a sweep-wide `g_nonblack_effects` counter: `smoke_one` bumps it on any non-zero frame sum and `run_effects_tests` asserts it ended positive after both the production and device roster passes — catches a total regression-to-black while still tolerating individually-black effects like RingShower.)*
41. ✅ `tests/test_transformers.h:136` — the ripple fast-reject thresholds are set to degenerate values that *disable* the reject path, so the rejection logic is never exercised and composition is not asserted to differ from a single entity. → Add a real-threshold case and a divergence assertion. *(Added `test_ripple_threshold_reject_path`: calls `prepare_thresholds()` and asserts an in-band peak point rotates while points past either bound (d&lt;d_min and d&gt;d_max) take both reject legs and return unchanged; and strengthened `test_transformer_spawn_applies_and_composes` to capture the single-entity output and assert the two-entity composition diverges from it.)*
42. ✅ `tests/test_styles.h`, `tests/test_generators.h` — `noise_warp`/`melt_warp` are exercised only with `noise=null` (identity path); the bound-noise production branch is untested, and `generators` has no multi-level-nesting stress for its reentrant scratch protocol. → Cover the bound-noise paths and a deep-nesting case. *(Added `test_noise_warp_bound_distorts` and `test_melt_warp_bound_noise_perturbs` — both bind a real `NoiseParams` via the production `sync_noise()` path and assert the warp takes its bound branch (displaces / diverges from the pure-drip result) while staying on the sphere; and `test_generate_deep_nesting_stacks_and_unwinds` drives 5 stacked `generate()` frames, asserting each enters strictly above the previous level's high-water, per-level sentinels survive, and the outermost scope rolls both arenas back.)*
43. ✅ `tests/test_reaction_graph.h:93` — `D_AVG` is pinned but nothing asserts `RD_N` equals the actual neighbor-table size, and the reciprocity tripwire (>50%) and CubemapLUT miss tolerance (¼) are loose for a near-exact lookup. → Add the self-consistency assert and tighten the thresholds. *(Added `test_table_shape_matches_constants`: `static_assert`s the declared `neighbors` shape equals RD_N×RD_K, plus a runtime check that the south-pole row holds real local neighbors — catching a generator that emitted fewer rows and silently zero-padded the tail. Tightened the reciprocity tripwire 50%→95% (measured 98.9%) and both CubemapLUT miss tolerances 25%→5% (measured 0 misses), adding an `[info]` line to the round-trip test so its counts are visible.)*
44. `tests/test_solids.h` — Catalan/Islamic solids are checked finite/non-empty but not asserted off-sphere, so a bug wrongly snapping `hankin`/`expand` vertices back to the unit sphere passes. → Assert some vertices are measurably off-sphere.

### Low

45. `core/static_circular_buffer.h:62` — the initializer-list overflow path calls `hs::log` unconditionally from a constructor that may run before `Serial.begin()` on device. → Gate behind `hs::debug` or document post-setup-only.
46. `core/static_circular_buffer.h:174-188` — internal `pop_*` paths decrement `count` with no `HS_CHECK`, relying on caller gating. → Add a symmetric guard or document the internal caller-guaranteed contract.
47. `core/concepts.h:296` — `ScalarFn = Fn<float(float),32>` is the lone callable-capacity choice with no rationale comment while every sibling is justified. → Add a one-line note on the 32-byte budget.
48. `core/inplace_function.h:139` — the static_assert message says "trivial POD/pointer closures" but the class accepts any nothrow-copyable callable. → Soften the wording.
49. `core/util.h:54` — `wrap(long,long)`/mixed-integer pairs route through `float fmod` and silently lose precision above 2²⁴. → Document that only the `(int,int)` overload is exact.
50. `core/easing.h:18-21` — cubic easings pass `t` unclamped and can return out-of-`[0,1]` for out-of-range input the sibling comments invite. → Confirm downstream clamps palette indices or note the easings are unclamped by design.
51. `core/waves.h:25` — `sin_wave` recomputes the loop-invariant `2π·phase` term per call. → Precompute in the capture if it ever moves to a hot loop.
52. `core/3dmath.h:23` — single-letter global `G` for the inverse golden ratio is easy to confuse. → Rename to `INV_PHI`.
53. `core/3dmath.h:1176` — Vector `slerp` antipodal fallback returns `v1` regardless of `t`. → Consider `t < 0.5 ? v1 : v2` so the degenerate path honors direction.
54. `core/geometry.h:304-314` — the float `y_to_phi` overload does not clamp `y` to row bounds (the int overload traps). → Document it is unchecked by design or add a debug assert for symmetry.
55. `core/geometry.h:802-810` — `get_antipode` negates two basis axes (preserving handedness) but the result no longer satisfies `make_basis`'s `w = cross(v,u)` relationship. → Verify no consumer recomputes `w` from `u,v`.
56. `core/spatial.h:25-119` — `AABB` is documented dead code ("no production consumer yet"). → Consider relocating to a `spatial_experimental.h` so the core header holds only wired-up code.
57. `core/spatial.h:240,261` — inconsistent `size_t` casting on `k` two lines apart. → Pick one.
58. `core/reaction_graph.h:113` — `lookup()` divides by the dominant-axis magnitude with no zero-guard, safe only because a unit `p` guarantees `≥1/√3`. → Note the unit-length precondition is load-bearing for the divide.
59. `core/color.h:1453` — `g_hue_seed` is a mutable static shared across all `GenerativePalette` instances (documented single-threaded; safe today because each WASM worker has its own linear memory). → Consider threading the seed through the constructor so the global is opt-in.
60. `core/color.h:1106` — calling the fixed coprime additive step (157 mod 256) a "low-discrepancy stride" overstates it; it is an R1 additive recurrence. → Reword.
61. `core/styles.h:166-173` — `sync_noise() const` mutates through bound state while sibling `sync_hue()` is non-const. → Add a note that constness reflects writes-through-bound-state vs own-members.
62. `core/memory.h:237` — `TriangularBitset::index()` uses a debug-only `assert` as the sole guard on its one memory-writing primitive; consistent with the hot-path doctrine but the highest-blast-radius `assert` here — on record as a deliberate device-side gap.
63. `core/effect_registry.h:82` — `reserve(64)` is a magic literal that silently under-reserves past 64 effects (harmless growth fallback). → Add a threshold comment (or `reserve(HS_EFFECT_COUNT)`).
64. `core/presets.h:83` — public mutable `entries`/`current_idx`/`prev_idx` force the runtime OOB trap in `get()`. → Make the indices private with `next()/prev()` the only mutators.
65. `core/canvas.h:39` — the `Effect` ctor runs two full `MAX_W*MAX_H` buffer fills despite `init()` existing to move heavy setup out of the ctor; defensible but worth noting.
66. `core/conway.h:343-349,360` — `relax` is listed once among buffered primitives and once as "no extra index buffers." → Clarify it is a no-scratch primitive in both lists.
67. `core/conway.h:489,497` — inner loops in `kis` reuse the name `i` (no shadow bug, outer loop closed) next to `offset + i`. → Rename to `k` as sibling operators do.
68. `core/mesh.h:539` — `classify_faces_by_topology` uses a 32-bit hash directly as the dense topology id with no tiebreak — safe for the fixed solid roster but the WASM mesh-editor is an untrusted boundary where a collision silently merges face classes. → Gate the finer keying behind the trusted/untrusted distinction.
69. `core/solids.h:1301` — `get_by_name` returns a live `cube()` after an "unreachable" `HS_CHECK(false)`; the "no silent wrong solid" guarantee depends on the macro being `[[noreturn]]`. → Use `__builtin_unreachable()`/a `[[noreturn]]` abort instead of a live fallback.
70. `core/hankin.h:214-219` — `dyn_idx` is captured *before* the matching `emplace_back`, reading as an off-by-one against the `size()-1`-after-push idiom elsewhere. → Add a "pre-push capture is deliberate" note.
71. `core/animation.h:1569` — comment cites "finding 419," violating the project's **no-finding-refs-in-comments** convention. → Drop the number and describe the behavior. (Convention fix.)
72. `core/animation.h:982-1003` — `Driver` re-reads `*speed_src_` every frame with no finite guard; a slider momentarily producing NaN permanently poisons the driven param via `wrap_t(NaN)`. → Guard with `std::isfinite`.
73. `core/filter.h:1182-1186` — `Screen::Trails::DecayPixel` array is named `ttls_` but holds position+ttl. → Rename to `points_` to match `World::Trails::items_`.
74. `core/transformers.h:36,137` — the template-level `params` member shares its name with `Entity::params`. → Rename to `template_params`/`default_params` to end the shadowing.
75. `core/transformers.h:288` — `ripple_transform` runs the full per-pixel `fast_acos`+`expf` even when amplitude is ~0 (between ripples). → Short-circuit on amplitude as `noise_transform` does.
76. `core/animation.h:2266` — the `Noise` 2²⁴-frame precision-freeze comment assumes `speed==1`. → Note the limit scales inversely with `speed`.
77. `core/filter.h:1614-1617` — `sample_bilinear_prev` hand-rolls a negative-x wrap instead of `fast_wrap` (correct, since `bx` can exceed the `fast_wrap` window). → Add a one-line note on why `fast_wrap` is inapplicable.
78. `hardware/pov_sync.h:1576` — `build_word_` is documented "written by ISR" but `seed()` also writes it once from the foreground at boot (benign, pre-attach). → Note the boot write.
79. `hardware/pov_segmented.h:659` vs `683` — `submitFrame()`'s `bool` return is checked in `render_black` but ignored in `render_column`. → Mark `submitFrame` `[[nodiscard]]` and `(void)`-cast the intentional discard.
80. `hardware/hd107s_frame.h:194` vs `147` — `packPixel` uses a stripped `assert` (hot ISR path) while `load()` uses always-on `HS_CHECK`. → Add a cross-reference so a future editor doesn't "harmonize" an always-on check onto the ISR hot path.
81. `effects/HopfFibration.h:129-131` — the stereographic pole comment implies the finite `1/eps` magnitude matters, but it is immediately discarded by `normalized_or`. → Clarify it only keeps pre-normalize direction well-defined.
82. `effects/Liquid2D.h:64,125` — `cycle_phase` uses a `Driver(wrap=false)` then a manual `fmodf`, two different wrap strategies in one init. → Let the Driver wrap it or comment why the manual wrap is preferred.
83. `effects/BZReactionDiffusion.h:24` — comment calls the lattice cells "Voronoi boundaries," falsely associating it with `Voronoi.h`. → Reword to "cell boundaries between lattice nodes."
84. `effects/Raymarch.h:107-194` — `shadeBlinnPhong` takes distinct `light_dir`/`view_dir` but every caller passes the same headlight vector. → Collapse to one `eye_dir` until a non-headlight caller exists, or keep per the doc.
85. `effects/Voronoi.h:159-164` — the border-color comment invokes premultiplied semantics that the final `color*alpha` multiply makes moot (any alpha-0 color yields black). → Simplify the comment.
86. `effects/IslamicStars.h:95` — the "every path must shuffle before draw" invariant on `palettes_slots` is enforced only by convention, unlike `HankinSolids`' guarded analog. → Add a guard or cross-reference.
87. `effects/MobiusGrid.h:93` — `phase_frame_` is advanced in lockstep with `timeline.step()` by assumption; an extra/missing step silently desyncs the phase. → Note the precondition or derive from a timeline-owned counter.
88. `effects/GnomonicStars.h:104` — the per-frame `fib_spiral` recompute is justified as "sim-targeted," but the effect is registered for device resolutions like every other. → Enforce the sim-only assumption or tighten the comment.
89. `effects/ShapeShifter.h:48-51` — **(carry-over, prior #9)** `Count` up to 128 dispatches ~256 SDF rasterizations/frame (both a Plot and a Scan pass per layer) with no per-frame budget cap (author-flagged). → Gate `Count` by resolution / clamp the effective loop against a pixel budget if it ships to a device.
90. `effects/Moire.h:62,169` — `density`'s in-class struct default (10.0f) is dead because `init()` overrides it per-resolution. → Note the override on the struct field.
91. `effects/Comets.h:258-262` — the class `@note` claims Comets applies `Screen::AntiAlias`, but it uses an empty `Pipeline<W,H>` and deliberately omits AA (the note appears copied from the `ChaoticStrings` sibling). → Correct the per-effect note.
92. `effects/Thrusters.h:104,188` — `t_global` is incremented at the end of `draw_frame`, so consumers see last frame's counter (documented invisible offset). → Increment at the top to remove the latent foot-gun.
93. `effects/Dynamo.h:144` — `color_wipe()` captures `&palette_boundaries.front()` whose validity depends on the `palettes.is_full()` (cap 16) guard never letting the cap-15 boundary buffer wrap onto a live slot. → Cross-reference that coupling at the capture site.
94. `effects/RingSpin.h:99` — all four rings spawn with the same `Y_AXIS` normal and diverge only over time, so the first frames render coincident (likely intentional, undocumented in an otherwise heavily-commented file). → Add a one-line note.
95. `effects/PetalFlow.h:220` — up to 64 rings × W/2 samples is the heaviest curve effect at 288×144; scratch handling is correct. → Profile the per-column budget on device.
96. `effects/MeshFeedback.h:151` — the `is_bound()` guard protects a state `init()` guarantees never occurs; harmless defensive guard, noted.
97. `effects/DreamBalls.h:247` — the `faces.size()/2 == E` closed-manifold assumption couples the preset table to the edge-count math and traps at `push_back` rather than a clear validation point (fail-fast by design); noted.
98. `targets/wasm/wasm.cpp:725` — `paramViews.reserve` is grouped in comments with the must-not-realloc view buffers, but it backs a freshly-built `val::array` (amortization only). → Clarify it is not part of the view-backing contract.
99. `targets/wasm/wasm.cpp:565` — `setParameter` returns the raw `updateParameter` bool so JS cannot distinguish "unknown name" from "rejected/clamped write." → Document the semantics like `setClip`/`setResolution`.
100. `targets/wasm/wasm.cpp:90` — the `tooling_scratch` "scratch outlives its op" contract is enforced only by WASM single-threading with no guard. → Add an assert if/when workers are introduced.
101. `targets/wasm/wasm.cpp:1010,1194` — `getMaxBounds` is gated behind `HS_WASM_DEV_BINDINGS` but no preset/Justfile defines it. → Document the build knob.
102. `CMakeLists.txt:44-72` — the debug 64 KB `STACK_SIZE` override wins only by link-line ordering; reordering the two `target_link_options` blocks would silently revert to 8 KB. → Set `STACK_SIZE` only in the build-type blocks to remove the order dependency.
103. `scripts/generate_luts.py:53` — `emit_array` emits a trailing comma and relies on an unpinned manual `clang-format` step to reach committed style. → Emit repo style directly or assert the format step ran.
104. `scripts/wasm_smoke.mjs:301` — `spline_catmull_rom_tangents` is checked only for finiteness while the cubic evaluators get endpoint pinning, so an argument transposition on this function passes. → Pin one known tangent value.
105. `scripts/capture_screenshots.mjs:119` — a resolution-resolve failure silently falls back to the app default and exits 0, so a gallery can be captured at the wrong resolution and pass. → Treat it as at least a loud summary warning.
106. `tools/teensy_gate.py:288` — `load_budgets` strips `//` `/* */` with naive regexes that would corrupt such sequences inside a future JSON string value (none today). → Use a JSONC-aware parser.
107. `tools/teensy_budgets.json:57` — the Phantasm `reaction_graph` FLASH invariant is `UNVERIFIED` against a real fitting ELF (the overflowing build is discarded). → Confirm once Phantasm fits, as the comment states.
108. `tools/teensy_warnings.py:36` — the warning ratchet keys on raw message text, so a compiler upgrade that re-words a warning reads as "new"; documented baseline-regen mitigates, noted.
109. `daydream/driver.js:316,471` — `==`/`!=` used for `e.key` and `stepFrames` comparisons where the codebase uses strict equality. → Switch to `===`/`!==`.
110. `daydream/driver.js:458` — the spiral-of-death catch-up clamp `0.25` is an inline magic number. → Extract a named constant with a rationale.
111. `daydream/driver.js:115` — `LABEL_VISIBILITY_COS` is described as a "cosine cutoff" but is a framing ratio exact only at the canonical camera distance. → Reword.
112. `daydream/daydream.js:131` — `refreshPixelView` relies on the load-bearing "non-detached ⇒ not stale" assumption (Emscripten growth detaches in place). → Add a one-line note.
113. `daydream/sidebar.js:100` — `setEffects()` rebuilds all button DOM via `innerHTML=''`, contradicting the class docstring's "reorder, not recreate" promise (which holds only for `sortBy`). → Clarify the docstring distinguishes the two.
114. `daydream/gui.js:74` vs `state.js:181` — duplicated `parseFloat(v.toFixed(4))` URL-serialization logic that can drift. → Extract a shared helper.
115. `daydream/segment_controller.js:567,599` — `updateStats` derives length from `this.timings` and `this.count` (kept in lockstep by `create()`). → Derive both from one source.
116. `daydream/recorder.js:201` — `_warnedNoRequestFrame` is used but never declared in the constructor, against the codebase's "all state visible" convention. → Initialize it.
117. `daydream/recorder.js:147` — `MediaRecorder` is passed `mimeType:''` when no codec is supported; some engines throw on an explicitly-empty option. → Omit the key to take the UA-default path by omission.
118. `daydream/segment_controller.js:471` / `segment_worker.js:162` — the row-wise blit loop and even-rounding helper are duplicated near-identically across the boundary. → Extract shared helpers into `segment_layout.js`.
119. `daydream/segment_worker.js:118` — `setResolution` guards on `=== false`; a binding returning `undefined` would pass through and update state for a size the engine may not have applied (holds today via the bool contract). → Assert a boolean or document the coupling.
120. `daydream/tools/solid_codegen.js:36` — `formatFloat` uses `val.toString()`, which emits scientific notation (`1e-7f`) for extreme magnitudes that does not compile as a C++ literal. → Route through `toFixed` like the sibling formatters.
121. `daydream/tools/spline_math.js:89` — `formatFloatCpp` strips trailing zeros before fixing the dangling dot, a fragile ordering. → Adopt the clearer two-step trim used in `lissajous_math.js`.
122. `daydream/tools/palette_math.js:328` — `GenerativePalette.get` interpolates in linear light; confirm the engine's `BakedPalette::get` interpolates in the same domain (linear vs sRGB) or soften the "matches the engine" claim. (Verify/document only — see the standing rejection in §5; do **not** re-interpret the sRGB-authored coefficients in OKLab.)
123. `daydream/tools/slider.js:67` — no validation that `min<max`/`step>0`; a misconfigured `scale=0` yields an inert slider with no error. → Assert the contract.
124. `daydream/tools/shared.js:98` — `initScene` reads `clientWidth/Height` with no null check on `getElementById`, throwing an opaque error on a wrong id (other modules guard this). → Add a null check.
125. `daydream/tests/*` — several test names/comments cite review-finding numbers ("finding 14", "finding 291"), violating the project's no-finding-refs convention. → Strip them from titles, keep in commit subjects. (Convention fix.)
126. `daydream/tests/color_parity_wasm.test.js`, `spline_math_wasm.test.js` — the parity guards run only if the WASM module is built; an absent/stale module yields a green run that silently skipped the comparison. → Ensure CI fails if the parity files are not collected/run.
127. `daydream/tests/segment_controller.test.js`, `segment_worker.test.js` — the WASM engine is mocked in both, so a drift between `FakeEngine`'s contract and the real `HolosphereEngine` method surface is not caught. → Add one thin real-engine integration test pinning the contract.

### Rejected (preserved from prior review — do not re-raise)

128. ❌ `core/color.h` — **Authoring VIBRANT `GenerativePalette` keys at the exact per-hue cusp introduced visible gradient banding and was reverted.** `get()` peaks chroma at mid-L *between* keys, so cusp-authored keys pushed that mid-segment peak past the gamut boundary; the in-gamut→chroma-clip switch then banded the gradient. The banding-safe chroma (~0.78·cusp) sits about where the fixed `kChromaPeak·sin(πL)` ceiling already is, so exact-cusp bought almost no vividness for the artifact. The companion idea — running `ProceduralPalette`'s cosines in OKLab — was also rejected: the 21 coefficient sets are sRGB-authored, so re-interpreting them as OKLab L/a/b drives every channel far out of range. (Standing rejection; finding 122 above is the *verify/document* descendant, not a re-proposal of this change.)

---

## 5. Reconciliation With the Prior Review

The thirteenth audit (graded A−, fifteen reviewers) reported the same overall profile this
pass reaches independently: no Critical/High defect, an A/A+ engine, and a gap-to-A that lives
at the test-coverage frontier. **The tree confirms its fix list is essentially fully
resolved:** the one genuine open functional defect that pass surfaced — the MeshFeedback
"Pause halts the shape carousel" guarantee being only half-implemented (its Finding 1, High) —
is closed (the re-armed hold timer now carries the `animationsPaused()` gate; this pass's
independent review of `MeshFeedback.h` found the pause path correct). The slerp `d==0`
boundary, the `update_hankin` stale-view asymmetry, the unclamped `MobiusFlow::progress`, the
borrowed-mode `transform` foot-gun, the unused `ParticleSystem` `max_delta` clamp, the
non-`explicit` `Quaternion`/`Color4` ctors, the `ChromaticShift` wrap divergence, the
`util::wrap` NaN guard, the bare-`volatile` segmented counters, the constructor-time hardware
I/O, the unguarded `#gui-container` append, the recorder mid-resize blit, and the `slider.js`
`innerHTML` XSS hole all carry verified-resolved dispositions in the current source.

**Carry-over open items** (independently re-confirmed open this pass, preserved with their
prior identity):

- Prior Finding 9 → **this pass Finding 89** — `ShapeShifter`'s unbudgeted per-frame
  rasterization cost (`Count` ≤ 128 → ~256 raster passes/frame, dual Plot+Scan), acknowledged
  in-code but still shipped.
- Prior Finding 10 → **this pass Finding 4** — 21–22 of 27–28 effects have no per-effect
  output invariant; smoke + determinism cannot catch a wrong-but-non-black regression.

**Standing rejection carried forward** (preserved so it is not re-litigated):

- Prior Finding 28 → **this pass Finding 128 (❌)** — the VIBRANT `GenerativePalette`
  exact-cusp authoring (gradient banding) and the `ProceduralPalette`-cosines-in-OKLab idea
  were both tried and rejected. This pass's Finding 122 is a *verify/document* descendant
  about the JS tool's interpolation domain, **not** a re-proposal of the rejected engine
  change.

**Where this pass adds information beyond the thirteenth audit:** a sharper accounting of the
JS browser-glue coverage gap (this pass elevates `driver.js`/`daydream.js`/`recorder.js`
lifecycle/`sidebar.js` being untested to High, Findings 1–3 — previously folded into the
general daydream-glue note), two latent Phantasm sync-timing items worth a confirming test
(the beacon mod-64 resync precondition, Finding 6; the same-tick burst/fold `j`-inference
ordering, Finding 7), the `SplineFlow` 240 KB device-arena allocation lacking a `static_assert`
(Finding 5), the unclamped `fib_spiral` `acosf` (Finding 13), the per-frame fade-multiply
down-bias (Finding 14), and two of the codebase's own no-finding-refs-in-comments convention
violations (Findings 71, 125). None of these contradict the prior audit; they refine the same
diagnosis.

This pass numbers its own fix list from 1 (the established per-pass convention); readers
tracking history should treat §4 as the current-tree superset and the thirteenth audit's `✅`
entries as closed.

---

## 6. Overall Appraisal — Magnitude of the Achievement

### What this project actually is

Holosphere is one C++ rendering engine that drives three radically different targets from a
single source of truth: a 600 MHz bare-metal Teensy 4.x painting a sphere with a strip of LEDs
spinning at 480 RPM; a four-board distributed sibling (Phantasm) that keeps four independent
CPUs phase-locked over a single wire with no shared clock; and a WebAssembly build that
reproduces the identical imagery in a browser as a design tool. On top of that substrate sit
27 original generative effects — reaction-diffusion on a Fibonacci lattice, Hopf fibrations,
spherical harmonics, Conway-operator polyhedra dressed in Hankin/Islamic star patterns,
Möbius-warped feedback, SDF ray-marching — all rendered in 16-bit linear light with OKLCH
perceptual palette interpolation, on a device with ~335 KB of usable RAM and no heap to spend.

### The magnitude of the achievement is large, along an unusual axis

The rare, hard thing here is not any single effect; it is the **engineering coherence** that
lets a firmware target, a distributed real-time system, and a browser application share one
codebase without the abstractions either leaking or collapsing. Eighteen independent reviewers
— most of whom saw only their slice of the tree — converged on the same verdict by the same
route: the way you find a flaw in this code is to locate the one place a self-imposed rule
*isn't* honored. That is the signature of a codebase whose discipline is real rather than
aspirational. Compile-time `<W,H>` parameterization gives the microcontroller zero-overhead
specialization while the simulator instantiates a different resolution from the same templates.
A partitioned arena with RAII scratch scopes replaces the heap entirely, every allocation seam
trapping and every budget documented. A fail-fast doctrine (`HS_CHECK`) is applied with genuine
discipline about *where* — cold seams trap, the per-pixel hot loop never does — and then
*verified* by a death harness that re-execs the binary and asserts each trap fires with the
exact illegal-instruction signal. That last move is the tell: most projects that talk about
fail-fast never prove their traps work; this one tests its own safety net, and proves the net
is load-bearing rather than decorative.

### Technical merit, relative to peers

- **Versus hobbyist / maker LED-art:** not in the same category — several tiers above. The
  typical artifact there is an Arduino sketch with `delay()`s, 8-bit sRGB blending, global
  mutable state, and no tests. Holosphere blends in linear light, eliminates the heap, templates
  the whole pipeline, mocks FastLED bit-faithfully so the host build is the device build's twin,
  and ships a sanitizer-instrumented native suite with a death harness and CI that drives every
  effect at both resolutions.
- **Versus professional offerings:** competitive on engine architecture and *ahead* on color
  correctness (linear-light + OKLCH + chroma-preserving gamut clipping, more correct than most
  shipping LED products) and testing rigor; the memory model, determinism contract, and
  ISR/concurrency reasoning are at or above commercial-firmware standards. Behind only where a
  solo author cannot match a vendor: fleet management, content-tooling ecosystems, and hardware
  productization.
- **Versus academic research:** the individual algorithms are largely known art (OKLab is
  Ottosson's, the SDF/sphere-tracing vocabulary Quílez's, reaction-diffusion/Conway/Hankin each
  have literatures). The originality is in the *integration* — running this entire menagerie in
  real time, deterministically, on a microcontroller, against a pixel-faithful browser twin,
  with a one-wire distributed-clock protocol holding four boards together and a
  death-harness-verified safety contract underneath. The Phantasm 1-wire sync alone is
  research-grade systems engineering and a credible systems-paper artifact.
- **Versus the state of the art:** for the specific niche of *synchronized multi-board
  full-color POV spheres rendered from one perceptually-correct templated engine with a
  bit-identical browser twin*, this appears to **be** the state of the art. There is no obvious
  public peer doing all of it at once.

### Artistic merit

The 27 effects are genuinely strong and *curated* rather than padded. The color pipeline
(perceptual interpolation, shortest-arc hue, gamut-preserving clipping, soft alpha in linear
light) is doing real aesthetic work that cheaper pipelines cannot, and it shows up exactly where
the README claims (soft gradients, multi-layer compositing). The medium itself — full-color
imagery floating on the surface of a physical spinning sphere — is striking, and the browser
simulator collapses an otherwise expensive iteration loop. This is the work of someone who is
simultaneously a serious systems engineer and a serious visual artist, and the two halves
reinforce rather than compromise each other.

### The honest gaps

The achievement is real and large; the gaps are narrow and known. They are: (1) the JS
browser-glue layer (`driver.js`, `daydream.js`, `recorder.js` lifecycle, `sidebar.js`) is
under-tested relative to the C++ core; (2) ~22/28 effects are behaviorally smoke-tested only, so
a visual regression can slip through green CI; (3) the device firmware budgets can't be verified
in automation, leaving a few large arena allocations confirmed only by reasoning; and (4) a
small set of latent timing/numeric edges (the two sync-protocol items, the `fib_spiral` clamp)
would each benefit from one confirming test. None undermines the work — they are the finishing
items on an already-exceptional engine.

### Bottom line

Holosphere is an **A− project whose distance to an A is almost entirely a coverage frontier,
not a quality deficit.** As an act of sustained, disciplined engineering across firmware,
distributed real-time systems, graphics, color science, a WASM bridge, and a web client — held
to a self-imposed standard that it then *verifies* rather than assumes — it is an unusually
accomplished piece of work. It comfortably exceeds its hobbyist peers, stands with professional
firmware on the axes that count, and in its sync subsystem and its sim/device-equivalence rigor
reaches genuinely state-of-the-art territory for the niche. The achievement is not a clever
effect; it is that, across three targets and two languages, none of the seams show.

---

*Review generated by independent multi-agent audit of the full codebase (README + ~18
component sub-reviews), then reconciled against the prior review (§5). No Critical or
High-severity correctness, memory-safety, or concurrency defect was substantiated. §4 is the
complete set of actionable items found.*
