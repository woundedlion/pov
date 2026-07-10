# Holosphere / daydream — Code Quality Review

**Scope.** A full two-repository audit of the Holosphere C++ rendering engine + firmware
and the daydream web simulator. Every in-scope source file was read by a dedicated
reviewer, graded against the architecture the README documents, and every candidate
finding was re-derived by an independent validator before inclusion. Out of scope by
request: `core/engine/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`,
`core/math/rotate.h`, and all vendored code (`core/vendor/`, `three.js/`,
`node_modules/`). The generated `holosphere_wasm.js`/`.wasm` artifacts in daydream are
build output, not reviewed.

**Method.** 21 component reviewers → 26 candidate findings → 26 independent per-finding
validators → **21 confirmed findings**. Five candidates were rejected or absorbed as
already-handled/by-design.

## Verdict

**Overall grade: A− (Exceptional).**

This is among the most disciplined codebases of its kind. Across ~85 hand-written source
files in two languages there is **not a single confirmed high or critical defect** — the
entire surviving finding set is 2 medium and 19 low, and most of those are documentation
drift, missing-coverage nits, or hardening niceties rather than live bugs. The engineering
philosophy stated in the README (16-bit linear light, compile-time `<W,H>` specialization,
arena-only allocation, fail-fast `HS_CHECK`, single-writer ISR discipline) is not just
described — it is *upheld consistently at every seam*, defended by compile-time
`static_assert`s, and verified by an unusually rigorous test suite (out-of-process death
tests, cross-run determinism, host-testable protocol cores, WASM runtime smoke, and
lockstep-drift-proof golden cross-checks).

### Grades by dimension

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness** | A | Math matches the documented contracts across projections, quaternions, CSG, physics, and protocol codecs; edge/pole/antipode/seam cases are handled and tested. The lone substantive miss is one arc-length register scale factor. |
| **Memory Safety** | A | Bump-arena model with wrap-proof subtractive bounds arithmetic, `static_assert` overflow ceilings, and debug-only use-after-free generation tracking. No render-path heap. One latent indeterminate-bytes-on-empty-copy nit. |
| **Concurrency & ISR Safety** | A | Rigorous single-writer ISR model on the Phantasm sync path (edge ISR = pure publisher, flywheel ISR = sole owner); the double buffer and mailbox handoff are race-argued and fused. One narrow worker-fault/retry-timer race in the simulator. |
| **Performance** | A | Exemplary hot-path discipline: split trig LUTs, branchless wraps, division-free rotate, measured-error `fast_*` approximations, zero-overhead template specialization. One exact-`atan2f` slip in a per-fragment shader. |
| **Architectural Elegance** | A | Clean domain layering (SDF/scan/canvas/shading; state/driver animation split; pure cores vs device shells; DOM-free logic vs glue). Purely-functional arena pipeline with polarity contracts documented at each seam. |
| **Portability** | A | Careful LP64-vs-32-bit wrap parity, `__FINITE_MATH_ONLY__` guards, bit-exact FastLED host mocks, and a native/WASM/device abstraction that keeps the simulator bit-identical to hardware. |
| **Error Handling** | A | Fail-fast/fail-dark applied deliberately: cold seams trap, hot paths use stripped asserts backed by cold traps, transient conditions get bounded handling. A couple of guard-ordering / boundary-completeness gaps. |
| **API Design / Interface Expressiveness** | A− | Intention-revealing interfaces (`normalized` vs `normalized_or`, named factories, compile-time `StaticPalette` composition, single-source param marshaling). Minor sibling-naming inconsistencies. |
| **Maintainability** | A− | Every non-obvious decision is documented at its site; naming is precise. A little template duplication and one convention-consistency gap. |
| **Documentation** | A− | Doc coverage is extraordinary (a 2,000-line README that reads as a spec, plus dense load-bearing doxygen). Docked for a handful of small code/doc drifts surfaced below. |
| **Testing & Verification** | A− | Death harness, determinism pass, meta-tests pinning the roster, host-tested protocol cores, WASM smoke, external-canonical goldens. Gaps are a few uncovered branches and one tautological assertion. |

---

## Prioritized fixes

All 21 confirmed findings, most-significant first. Numbering is sequential across the whole
list; each is annotated with `[repo] file:line` and severity. There are no P0 (critical/high)
items.

### P1 — Correctness & performance (medium)

1. ✅ **[Holosphere] `core/render/plot.h:1130` — Ring `v1` arc-length register uses the wrong latitude scale (medium, correctness).** `Ring::sample` sets `arc_scale = sinf(work_radius)` and writes `v1 = theta * arc_scale`, but the ring's colatitude is `theta_eq = work_radius*(PI/2)`, so the true per-radian arc length is `sinf(theta_eq)` — already computed as `r_val`. At the equator Ring reports a perimeter `2π·sin(1)` (~16% short) and disagrees with `DistortedRing` (which accumulates true geodesic length) by up to ~48% for the identical ring, so a shader keying a gradient off the documented "arc length (radians)" register renders at a different density on `Plot::Ring` than on sibling primitives. Fix: use `r_val` in both `Ring::sample` overloads (lines 1130, 1197) and the close vertices (1159, 1229), and correct the local doc at line 1116.

2. ✅ **[Holosphere] `effects/Raymarch.h:144` — exact `atan2f` on the per-fragment ray-march hot path (medium, performance).** `frag_fn` runs once per lit torus pixel across 26 tori and computes the palette-lookup azimuth with libm `atan2f(loc.z, loc.x)` (~80–120 cycles vs ~8 for `fast_atan2`). It is the only transcendental in the fragment path and drives only a hue LUT index, which tolerates `fast_atan2`'s ~0.1° error — and every other azimuth-for-coloring site in the engine (`sdf.h`, `SphericalHarmonics.h`) already uses `fast_atan2`. Fix: `fast_atan2(loc.z, loc.x)`; it returns `[-π,π]` so the surrounding remap is unchanged.

### P2 — Correctness, safety & robustness (low)

3. ✅ **[Holosphere] `core/animation/timeline.h:193` — `add_get` pin-check comment asserts a false invariant, deferring a trap.** The loop admits any repeating/infinite predecessor because "repeating or infinite predecessors are never removed"; but a repeating animation that is later `cancel()`ed *is* removed, and the ensuing compaction relocates the pinned successor into `move_into`'s `HS_CHECK(!handled)` trap. Outcome is a fail-fast trap (no corruption), so low — but the justifying invariant is inaccurate. Fix: correct the comment to note the `cancel()` path and that the real backstop is the `move_into` check.

4. ✅ **[daydream] `gui.js:235` — step-snap uses exact `!==` against the raw value, so float noise marks an on-grid value as clamped and fires a redundant URL rewrite.** After snapping to a step multiple, `valClamped = val !== raw` becomes true for an already-on-grid input (e.g. `0.3` → `0.30000000000000004`), scheduling an unnecessary `history.replaceState` on load. Live for real stepped controls (e.g. `columnFillOverlap`, step 0.01). Fix: compare with a step-scaled tolerance, or only rewrite when the value was truly clamped/rejected.

5. ✅ **[daydream] `segment_controller.js:513` — `onWorkerFault` does not cancel a pending boot-retry timer, so a latched fault can be silently rebuilt away.** If a transient bare-`Event` boot failure arms the 250 ms retry and a deterministic fault (version mismatch / `initFailed` / worker throw) latches before it fires, the still-armed timer calls `create()`→`destroy()` and clears the "no auto-restart by design" latch. Narrow window, self-limiting, cosmetic (overlay flicker). Fix: add `this.clearRetryTimer();` alongside the other `clear*` calls at the top of `onWorkerFault()`.

6. ✅ **[Holosphere] `core/engine/inplace_function.h:158` — copy/move-from-empty leaves `storage_` indeterminate, contradicting the file's stated invariant.** The default/nullptr ctors value-initialize `storage_{}` "so a future op that reads an empty function's buffer is well-defined," but the copy (158) and move (161) ctors don't, and the empty-vtable ops are no-ops — so `auto b = a;` on an empty function yields indeterminate bytes. Harmless today (no op reads an empty buffer). Fix: add `storage_{}` to both ctors, or narrow the comment to the default/nullptr ctors.

7. ✅ **[Holosphere] `core/math/3dmath.h:1215` — `fast_expf` saturation guard runs after a potentially-UB float→int cast, and its documented `x<=0` domain has no debug assert.** `int i = (int)fi;` precedes `if (i < -126)`, so an astronomically negative `x` (unreachable with current bounded callers) would hit signed-overflow UB before the guard; and unlike sibling hot-path helpers (`operator/` carries `assert(s != 0)`), the documented `x<=0` precondition is unenforced. Fix: move the range check into float space (`if (fi < -126.0f) return 0.0f;`) before the cast, and add `assert(x <= 0.0f)`.

8. **[Holosphere] `targets/wasm/wasm.cpp:1337` — `spline_cubic_slerp` aborts the whole WASM module on finite-but-degenerate control points.** `eval_cubic_spline` gates only on `all_finite()`, but a finite zero-length control point drives `slerp → normalized()` into its always-on `HS_CHECK`, trapping the module — exactly the JS-boundary abort the `all_finite()` docstring claims to prevent, and inconsistent with the sibling `spline_cubic_fast` (which degrades via `normalized_or`). Not reachable from `splines.html`'s unit-sphere inputs; a defensive-boundary/doc-accuracy gap. Fix: reject degenerate control points at the boundary (or route the slerp blends through `normalized_or`), and correct the docstring.

9. **[daydream] `tools/mobius.html:647` — fine keyboard arrow step (0.05) is swallowed by `snapComplex`'s 0.1 zero-band, pinning zeroed Möbius coefficients.** Un-shifted arrows step 0.05, then `snapComplex` zeros any `|v| < 0.1`; starting from a defaulted-zero coefficient, each fine press recomputes `0.05` and snaps back to 0, so only Shift+Arrow (0.2) can move it — defeating the keyboard-accessibility path the handler exists to provide. Fix: raise the arrow step above the band, pass a smaller threshold for keyboard input, or skip `snapComplex` on discrete key events.

### P3 — Test coverage, maintainability, API & documentation (low)

10. **[Holosphere] `tests/test_plot_scan.h:666` — the Ring `v1` arc-length test is tautological.** It recomputes `arc_scale = sinf(wr)` and asserts `v1 == theta*arc_scale`, checking the code against a copy of its own (wrong) formula — so it passes despite finding #1. Fix: assert Ring's full-perimeter `v1` against the true great-circle arc (`2π·sinf(theta_eq)`) or against `DistortedRing::sample`'s `v1`.

11. **[Holosphere] `effects/Liquid2D.h:189` — the hand-derived glitch-lens rational map has no unit-norm/pole regression test.** `apply_glitch_lens` is a degree-3 rational map (hemisphere fold + triple-θ + `R²<1e-6` pole branch) on the live per-pixel path; it is currently correct but a sign/coefficient slip would be invisible to the positive-frame-sum smoke harness, and the existing white-box seam doesn't reach it. Fix: add a unit test asserting `|apply_glitch_lens(v)| ≈ 1` over a spread of directions plus the pole-branch return.

12. **[daydream] `tests/spline_math.test.js:184` — `randomPointOnSphere`'s rejection loop is never exercised.** The sole test feeds values chosen to "avoid the `s>=1`/`s===0` reject," so an inverted reject condition (accepting `s>=1` → non-unit point, or looping on `s<1` → hang) would pass. Fix: add a sequence whose first pair yields `s>=1` (rejected) then an accepted pair, asserting the result is unit-length.

13. **[Holosphere] `core/animation/trails.h:67` — `OrientationTrail` and `VectorTrail` are near-identical duplicated templates.** Both wrap `StaticCircularBuffer` with identical `record/length/get/clear/expire` forwarders, differing only in element type; the two `get()` docs have already drifted. Fix: unify as `template<typename T,int CAP> class Trail` with `using` aliases, preserving `VectorTrail::CAPACITY`.

14. **[Holosphere] `effects/FlowField.h:45` — relies on the default arena split with no compile-time footprint guard.** Its 600-particle pool is carved from the default `persistent_arena` without the `FOOTPRINT_BYTES` `static_assert` that sibling effects (HopfFibration, GnomonicStars) use to pin the same pattern — so a future retune of the pool or the default partition moves the failure from a build error to a runtime device trap. Fix: add the `FOOTPRINT_BYTES ≤ persistent-budget` `static_assert` plus a one-line note.

15. **[Holosphere] `core/engine/memory.h:801` — inconsistent empty-check naming: `ArenaVector::is_empty()` vs `ArenaSpan::empty()`.** The owning and borrow types in the same header use different spellings (`StaticCircularBuffer` also uses `is_empty()`), defeating generic code. Purely cosmetic (swaps produce loud compile errors, and `empty()` is defensible STL-view conformance). Fix: standardize on `is_empty()`, optionally keeping `empty()` as an alias on `ArenaSpan`.

16. **[daydream] `segment_layout.js:44` — error messages misattribute their source as `segment_worker:`.** All five `computeSegmentRange` throws prefix `segment_worker:`, but this pure-math module was deliberately factored out for Node testing and also runs in the controller's composite path with no worker on the stack, so a debugger is pointed at the wrong file. Fix: change the prefix to `segment_layout:` (or drop it) in all five throws.

17. **[Holosphere] `README.md:615` — the Fragment comment inverts `Color4()`'s default alpha.** It says "`Color4()`'s default is opaque," but `color.h:229` defines `Color4()` with `alpha(0.0f)` (transparent) — the value `shading.h` relies on to cull a fragment via `return Color4()`. The parenthetical even contradicts the correct clause before it. Fix: state that `Color4()` defaults to transparent black (alpha 0.0), or drop the parenthetical.

18. **[Holosphere] `README.md:1261` — hardware-ID decode formula omits the active-low inversion the code applies.** README §1 / §7.10 print `raw & (N-1)`, but `pov_segmented.h:417` is `(~raw0) & (N-1)` because the straps are `INPUT_PULLUP` active-low; taken literally the README formula yields ID 3 for the all-floating master, contradicting its own "all-floating = ID 0" note. Fix: print `(~raw) & (N-1)` in both locations, or define `raw` as the post-inversion logical value.

19. **[daydream] `driver.js:612` — the PiP camera copies the main camera exactly, contradicting its documented "front and back simultaneously" purpose.** `renderPip()` copies both `position` and `quaternion` each frame (differing only in square aspect), and the cull uniform uses the main camera position, so the PiP can never show the back hemisphere — the second render pass is wasted work. Fix: give `pipCamera` a genuinely fixed/180°-offset orientation, or delete the redundant path and update README §10.3.

20. **[Holosphere] `core/color/color.h:344` — `srgb_to_linear_interp` `@param` claims callers must clamp, but the body clamps internally.** The doc says "callers must clamp before calling"; the body unconditionally does `hs::clamp(...)` (required for float→int UB safety). A maintainer trusting the contract could remove the internal clamp and reintroduce the exact UB the inline comment warns against. Fix: reword the `@param` to state that out-of-range/NaN inputs are clamped.

21. **[daydream] `tests/mobius_transforms.test.js:179` — incorrect magnitude in a `cdiv` guard comment.** The comment reads `|q|² = 4e-7`, but `q = {re: 4e-4}` gives `1.6e-7`. The assertion still holds (the test is valid); only the stated value is wrong and would mislead anyone recomputing the threshold margin. Fix: change `4e-7` to `1.6e-7`.

---

## What the codebase does exceptionally well

- **A stated philosophy that is actually enforced.** 16-bit linear light reaches unbroken from canvas to SPI wire; `<W,H>` specialization is genuinely zero-overhead; arena allocation is explicit at every call site with documented ping-pong polarity; `HS_CHECK` guards cold seams while hot paths use compiled-out asserts backed by cold traps. These are not aspirations in comments — they hold at every seam, and the death harness proves the traps fire.
- **Compile-time invariants over runtime hope.** CSG span budgets, arena footprints, pipeline filter ordering, param-marshal roster order, preset field layout, and integer overflow ceilings are all pinned by `static_assert`, converting whole classes of would-be runtime failures into build errors.
- **Hardware correctness reasoned like a datasheet.** The Phantasm one-wire sync protocol (flywheel timebase, odd-only distance-2 symbol alphabet, epoch-counted playlist, fail-dark identity) is a genuinely sophisticated distributed-timing design, and its entire load-bearing core is factored into host-testable pure code and unit-tested without a Teensy.
- **Test engineering that resists rot.** Out-of-process death tests with a trap-shape sentinel, a determinism pass that dirties every reset global so each reset line is self-testing, meta-tests that pin the module roster bidirectionally, external-canonical color goldens (Ottosson OKLab), and lockstep-drift-proof cross-checks between JS ports and the C++ source.
- **Honest, load-bearing documentation.** The README is effectively a design spec, and the doxygen explains the *why* of non-obvious invariants (cache clean-not-invalidate, straight-alpha vs premultiply, view detachment under `ALLOW_MEMORY_GROWTH`, coprime hue cursors). The doc findings above are drift in an unusually thorough corpus, not gaps.

## Residual risk

Negligible for a project of this scope. The two medium findings are a register-semantics
inaccuracy (visual-density only, not a crash) and a per-pixel micro-optimization; the
remaining nineteen are latent hardening, missing-coverage, or documentation items with no
live failure in normal operation. Nothing here blocks deployment; the list is a polish
backlog, not a bug queue.
